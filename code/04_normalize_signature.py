#!/usr/bin/env python3
"""
HISTORICAL IMPLEMENTATION: not part of the audited default pipeline.
CPM coefficients are not validated cell-number fractions. See docs/METHODS.md.

04_normalize_signature.py -- put the cell-type signature S on a CPM scale.

WHY THIS EXISTS (this is not cosmetic)
--------------------------------------
build_signature.py writes S as the MEAN COUNTS PER CELL of each cell type:

    S[k, g] = (sum of counts of gene g over cells of type k) / n_cells_k

The four cell types have wildly different RNA content / sequencing depth.
Measured total counts in the current signature:

    Postmeiotic    100,091,477
    Meiotic         92,981,459
    Spermatogonia    9,236,448
    Somatic          6,955,472

Postmeiotic carries ~14x the counts of Somatic. If we feed this straight into
NNLS, the solver is free to explain the data by inflating the Postmeiotic
coefficient (its columns are literally ~14x larger), so the estimated
composition is biased away from the truth. Independent simulation from the
The deconvolution spec (docs/03-null-model.md, section 2f) quantified the damage:
CPM-normalising S improves proportion-recovery MAE by ~4.6x
(0.151 -> 0.033), and NOT normalising overestimates the Somatic fraction
by up to 0.26 in absolute proportion -- which is the same order as the
composition shift we are trying to measure, i.e. it would manufacture the
very signal the project is testing for.

WHAT IT DOES
------------
Rescales each row (cell type) to counts-per-million, so that the K profiles
are comparable in scale and the NNLS coefficients are interpretable as
"millions of transcripts contributed by cell type k". Nothing else changes:
gene order and column names are preserved exactly, and the output is written
next to the input as `<stem>_cpm.csv.gz`.

Usage
-----
    .venv/bin/python code/04_normalize_signature.py \
        --in  data/ref/signature_GSE112393.csv.gz \
        --out data/ref/signature_GSE112393_cpm.csv.gz
"""

import argparse
import sys

import numpy as np
import pandas as pd

TARGET_SUM = 1e6  # CPM


def main() -> int:
    ap = argparse.ArgumentParser(description="CPM-normalise a cell-type signature matrix")
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--target", type=float, default=TARGET_SUM,
                    help="per-cell-type row sum after scaling (default 1e6 = CPM)")
    a = ap.parse_args()

    S = pd.read_csv(a.inp, index_col=0)
    print(f"read {a.inp}: {S.shape[0]} cell types x {S.shape[1]} genes")

    if (S.to_numpy() < 0).any():
        print("ERROR: signature contains negative values; refusing to scale", file=sys.stderr)
        return 1

    raw_sums = S.sum(axis=1)
    print("\nrow sums BEFORE (raw mean counts per cell):")
    for k, v in raw_sums.items():
        print(f"  {k:14s} {v:>16,.0f}")

    dead = raw_sums[raw_sums <= 0].index.tolist()
    if dead:
        print(f"ERROR: cell types with zero total counts: {dead}", file=sys.stderr)
        return 1

    n_cells = S.shape[1]
    Sc = S.div(raw_sums, axis=0) * a.target

    print(f"\nrow sums AFTER (scaled to {a.target:,.0f}):")
    for k, v in Sc.sum(axis=1).items():
        print(f"  {k:14s} {v:>16,.1f}")

    if not np.allclose(Sc.sum(axis=1).to_numpy(), a.target, rtol=1e-9):
        print("ERROR: row sums are not exactly the target", file=sys.stderr)
        return 1
    if Sc.shape != S.shape or not (Sc.columns.tolist() == S.columns.tolist()):
        print("ERROR: gene space changed during scaling", file=sys.stderr)
        return 1

    Sc.to_csv(a.out, compression="gzip", float_format="%.6g")
    print(f"\nwrote {a.out}  shape={Sc.shape}")
    print("note: gene order and column names are byte-identical to the input")
    return 0


if __name__ == "__main__":
    sys.exit(main())
