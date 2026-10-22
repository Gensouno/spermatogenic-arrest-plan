#!/usr/bin/env python3
"""Stream atlas counts into source-batch x cell-class pseudobulks.

Source batches are not assumed to be independent animals. Keep all known
classes and record excluded Unknown cells. No dense cell-by-gene allocation.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GROUPS = {1:'Somatic', 2:'Somatic', 3:'Somatic', 4:'Somatic', 5:'Somatic',
          6:'Somatic', 8:'Spermatogonia', 9:'Meiotic', 10:'Postmeiotic', 11:'Postmeiotic'}


def main():
    started = time.monotonic()
    base = ROOT / 'data/raw/ref'
    dge = base / 'GSE112393_DGE.txt.gz'
    annotation = base / 'GSE112393_PerCellAttributes.txt.gz'
    out = ROOT / 'data/ref/batches'
    out.mkdir(parents=True, exist_ok=True)
    ann = pd.read_csv(annotation, sep='\t', skiprows=3)
    ann['class'] = ann.CellType.map(GROUPS)
    if not set(ann.CellType).issubset(set(GROUPS) | {7}):
        raise ValueError('Unknown annotation code')
    keys = sorted(set(zip(ann.loc[ann['class'].notna(), 'Dataset'],
                          ann.loc[ann['class'].notna(), 'class'])))
    lookup = {k:i for i,k in enumerate(keys)}
    assignment = np.array([lookup.get(k, len(keys)) for k in zip(ann.Dataset, ann['class'])])
    genes, sums = [], []
    with gzip.open(dge, 'rt') as stream:
        if stream.readline().strip().split('\t') != ann['#CellBarcode'].tolist():
            raise ValueError('Barcode ordering mismatch')
        for n, line in enumerate(stream, 1):
            gene, sep, values = line.rstrip().partition('\t')
            counts = np.fromstring(values, sep='\t', dtype=np.float64)
            if not sep or len(counts) != len(ann) or not np.isfinite(counts).all() or counts.min() < 0:
                raise ValueError(f'Invalid DGE row {n}')
            genes.append(gene)
            sums.append(np.bincount(assignment, weights=counts, minlength=len(keys)+1)[:len(keys)])
            if n % 5000 == 0:
                print(f'{n} genes; {time.monotonic()-started:.1f}s', flush=True)
    raw = np.stack(sums, axis=1)
    if len(set(genes)) != len(genes):
        raise ValueError('Duplicate reference genes require explicit aggregation')
    cells = np.bincount(assignment, minlength=len(keys)+1)[:len(keys)]
    meta = pd.DataFrame(keys, columns=['source_batch', 'cell_class'])
    meta['cells'] = cells
    meta['total_umi'] = raw.sum(axis=1).astype(np.int64)
    meta.to_csv(out / 'batch_classes.csv', index=False)
    np.savez_compressed(out / 'counts.npz', sums=raw, genes=np.array(genes), cells=cells)
    provenance = {'cells_total':len(ann), 'cells_unknown_excluded':int(ann['class'].isna().sum()),
                  'genes':len(genes), 'source_batches':int(ann.Dataset.nunique()),
                  'elapsed_seconds':round(time.monotonic()-started, 2),
                  'note':'Source batches are not verified independent donors. Enriched batches retained for sensitivity; no abundance inference from captured cell counts.',
                  'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (dge, annotation)}}
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(meta.to_string(index=False), flush=True)
    print(provenance, flush=True)


if __name__ == '__main__':
    main()
