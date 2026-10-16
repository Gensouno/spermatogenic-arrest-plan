#!/usr/bin/env python3
"""Build the cell-type signature matrix S from the GSE112393 single-cell atlas.

Why this is written the way it is
---------------------------------
The GSE112393 merged DGE matrix is 34,633 cells x ~20k genes. Loaded densely as
float64 that is 34633 * 20000 * 8 = 5.5 GB, and pandas would roughly double it
during parsing -> ~11 GB in a 16 GB machine. That is the single step most likely
to kill the pipeline.

So we never materialise the full matrix. We stream the file in chunks of *genes*
(rows), and for each chunk immediately collapse the 34,633 cell columns down to
11 cell-type sums using the PerCellAttributes annotation. Peak memory is then
chunk_genes * 34633 * 8 bytes (~110 MB at chunk=400) plus the tiny accumulator
(11 x chunk_genes).

Output: a K x G matrix of pseudobulk-mean expression per cell type, which is
exactly the S in  X_bulk ~ S . P  used for NNLS deconvolution.

Usage:
    python3 scripts/build_signature.py \
        --dge data/raw/ref/GSE112393_DGE.txt.gz \
        --anno data/raw/ref/GSE112393_PerCellAttributes.txt.gz \
        --out data/ref/signature_GSE112393.csv.gz \
        --collapse data/ref/cluster_collapse.csv
"""
import argparse
import gzip
import sys
import time

import numpy as np
import pandas as pd

# The 11 annotated cell classes in GSE112393 (from the PerCellAttributes header).
CELLTYPE = {
    1: "InnateLymph", 2: "Macrophage", 3: "Endothelial", 4: "Myoid", 5: "Leydig",
    6: "Sertoli", 7: "Unknown", 8: "SPG", 9: "Scytes", 10: "STids", 11: "Elongating",
}

# Default collapse: bulk deconvolution cannot reliably separate the four
# consecutive germ-cell states from each other, and the scientific question is
# about the germ/somatic composition shift, not about intra-germ resolution.
# Set --collapse none to keep all 11 classes.
DEFAULT_COLLAPSE = {
    "Somatic": ["InnateLymph", "Macrophage", "Endothelial", "Myoid", "Leydig", "Sertoli"],
    "Spermatogonia": ["SPG"],
    "Meiotic": ["Scytes"],
    "Postmeiotic": ["STids", "Elongating"],
    # Unknown is dropped: it carries no interpretable share and would absorb
    # residual variance without a biological meaning.
}


def load_annotation(path: str) -> pd.DataFrame:
    """Parse the PerCellAttributes file.

    Layout (verified against the real file):
      line 0: '###### Per-cell Attributes for All 34,633 Cells ...'
      line 1: '### Cell type: 1-InnateLymph; 2-Macrophage; ...'
      line 2: '### Each row is one cell'
      line 3: '#CellBarcode\tDataset\tCellType\tnGene\t...'   <- real header (27 fields)
      line 4+: data
    We only want the first three columns, so we name them ourselves instead of
    trusting the 27-field header.
    """
    df = pd.read_csv(path, sep="\t", skiprows=4, header=None,
                     names=["barcode", "dataset", "celltype"],
                     usecols=[0, 1, 2], dtype={0: str, 1: str, 2: int},
                     compression="gzip")
    df["celltype_name"] = df["celltype"].map(CELLTYPE)
    if df["celltype_name"].isna().any():
        bad = sorted(df.loc[df["celltype_name"].isna(), "celltype"].unique())
        raise ValueError(f"unmapped CellType code(s): {bad}")
    return df


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dge", required=True)
    ap.add_argument("--anno", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--group-out", required=True)
    ap.add_argument("--chunk", type=int, default=400, help="genes per chunk")
    ap.add_argument("--collapse", default="default",
                    choices=["default", "none"],
                    help="'default' merges the 11 classes into 4 bulk-resolvable groups")
    a = ap.parse_args()

    t0 = time.time()
    print(f"[1/4] reading annotation from {a.anno}", flush=True)
    anno = load_annotation(a.anno)
    print(f"      {len(anno):,} cells; celltypes: "
          f"{anno['celltype_name'].value_counts().to_dict()}")

    # Map each cell column to a signature row index.
    groups = anno["celltype_name"].tolist()
    if a.collapse == "default":
        lookup = {}
        for new, members in DEFAULT_COLLAPSE.items():
            for m in members:
                lookup[m] = new
        groups = [lookup.get(g) for g in groups]
    # Drop cells whose class we do not model (Unknown).
    keep = np.array([g is not None for g in groups])
    group_names = sorted({g for g in groups if g is not None})
    print(f"      modelling {len(group_names)} signature rows: {group_names}")
    print(f"      dropping {int((~keep).sum()):,} cells (no modelled class)")

    gidx = {g: i for i, g in enumerate(group_names)}
    col_group = np.array([gidx[g] if g is not None else -1 for g in groups], dtype=np.int32)

    # The DGE columns are matched to annotation rows BY POSITION, so the two
    # files must be in the same barcode order. Verified identical on the real
    # data; assert it here so a future re-download cannot silently mis-align
    # the signature (which would produce a plausible-looking but wrong S).
    with gzip.open(a.dge, "rt") as fh:
        dge_barcodes = fh.readline().rstrip("\n").split("\t")
    if dge_barcodes != anno["barcode"].tolist():
        raise ValueError(
            "DGE column barcodes are not in the same order as the annotation "
            "rows; the positional join would be wrong. Re-align explicitly.")

    # Stream the DGE gene-chunk by gene-chunk using pandas' native chunked reader.
    # Peak memory = chunk_genes * 34,633 * 4 bytes (float32) ~ 55 MB at chunk=400,
    # plus the accumulator K x G which is only a few MB. The one-hot (n x K) is
    # rebuilt once outside the loop.
    print(f"[2/4] streaming DGE in chunks of {a.chunk} genes ...", flush=True)
    K = len(group_names)
    ncell = len(col_group)
    oh = np.zeros((ncell, K), dtype=np.float32)
    valid = col_group >= 0
    oh[np.arange(ncell)[valid], col_group[valid]] = 1.0
    oh[~keep] = 0.0
    acc_ncell = oh.sum(axis=0).astype(np.float64)   # constant per chunk

    acc_sum = None          # K x G accumulated raw counts per cell type
    gene_names: list[str] = []
    nrows = 0
    header_n = None
    with gzip.open(a.dge, "rt") as fh:
        # The DGE header line lists ONLY the 34,633 cell barcodes (no leading
        # field for the gene column), while every data row has 34,634 fields
        # (gene name + 34,633 counts). Verified against the real file. So we
        # must NOT let pandas parse line 1 as a header -- skip it and name the
        # columns ourselves, taking field 0 as the index.
        # dtype must be given as a per-column dict: passing the scalar
        # np.float32 makes pandas attempt the conversion BEFORE index_col
        # extraction and it dies on the gene-name column with
        #   ValueError: could not convert string to float: '0610005C13Rik'
        col_dtype = {i: np.float32 for i in range(1, ncell + 1)}
        reader = pd.read_csv(fh, sep="\t", header=None, skiprows=1,
                             index_col=0, chunksize=a.chunk, dtype=col_dtype)
        for buf in reader:
            if header_n is None:
                header_n = buf.shape[1] + 1   # +1 for the index column
                print(f"      header lists {buf.shape[1]:,} cell columns")
                if buf.shape[1] != ncell:
                    raise ValueError(
                        f"DGE has {buf.shape[1]} cell columns but annotation "
                        f"has {ncell}; refusing to proceed with a misaligned join")
            nrows += buf.shape[0]
            gene_names.extend(map(str, buf.index))
            mat = np.ascontiguousarray(buf.to_numpy(dtype=np.float32))
            block = mat @ oh                   # chunk_genes x K
            add = block.T.astype(np.float64)   # K x chunk_genes
            if acc_sum is None:
                acc_sum = add
            else:
                # last chunk may be shorter than the rest
                acc_sum = np.hstack([acc_sum, add])
            if nrows % (a.chunk * 5) == 0:
                print(f"      {nrows:,} genes  ({time.time()-t0:.0f}s)", flush=True)

    if acc_sum is None:
        print("ERROR: no data read", file=sys.stderr)
        return 2
    print(f"[3/4] {nrows:,} genes x {K} cell types "

          f"({time.time()-t0:.0f}s)", flush=True)

    # Per cell type we keep the mean count per cell (a pseudobulk profile).
    mean = acc_sum / np.maximum(acc_ncell, 1)[:, None]
    out = pd.DataFrame(mean, index=group_names, columns=gene_names)
    out.index.name = "celltype"

    import os
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    out.to_csv(a.out, compression="gzip")
    print(f"[4/4] wrote {a.out}  shape={out.shape}")

    det = pd.DataFrame({
        "celltype": group_names,
        "n_cells": acc_ncell.astype(int),
        "detected_genes": (mean > 0).sum(axis=1),
        "total_counts": acc_sum.sum(axis=1).astype(np.int64),
    })
    det.to_csv(a.group_out, index=False)
    print(det.to_string(index=False))
    print(f"done in {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
