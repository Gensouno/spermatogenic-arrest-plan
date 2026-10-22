#!/usr/bin/env python3
"""Match atlas source batches to official GEO sample metadata."""
import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def main():
    path=ROOT/'data/raw/metadata/GSE112393_family.soft.gz'
    rows, current = [], None
    for line in gzip.decompress(path.read_bytes()).decode().splitlines():
        if line.startswith('^SAMPLE = '):
            current={'gsm':line.split(' = ',1)[1]}
            rows.append(current)
        elif line.startswith('^'):
            current=None
        elif current is not None and line.startswith('!Sample_title = '):
            current['source_batch']=line.split(' = ',1)[1]
        elif current is not None and line.startswith('!Sample_characteristics_ch1 = '):
            key, val=line.split(' = ',1)[1].split(': ',1)
            current[key]=val
    meta=pd.DataFrame(rows)
    ann=pd.read_csv(ROOT/'data/raw/ref/GSE112393_PerCellAttributes.txt.gz',sep='\t',skiprows=3)
    n=ann.Dataset.value_counts().rename_axis('source_batch').reset_index(name='captured_cells')
    merged=n.merge(meta,on='source_batch',validate='one_to_one',how='left')
    unresolved=merged.loc[merged.gsm.isna(),'source_batch'].tolist()
    if set(unresolved)-{'ST7','ST8'}:
        raise ValueError(f'Unexpected missing official metadata: {unresolved}')
    # Annotation ST7/ST8 do not match GEO titles 1nDepleted1/2. Do not invent
    # a one-to-one assignment without an authoritative alias table.
    meta[~meta.source_batch.isin(n.source_batch)].to_csv(
        ROOT/'results/reference_calibration/unmatched_geo_sources.csv',index=False)
    merged['included_ST_only_reference']=merged.source_batch.str.match(r'^ST[1-6]$')
    assert (merged.loc[merged.included_ST_only_reference,'cell type enrichment']=='none').all()
    out=ROOT/'results/reference_calibration'
    merged.to_csv(out/'official_atlas_sources.csv',index=False)
    # The original annotation header is the authority for class code meanings.
    with gzip.open(ROOT/'data/raw/ref/GSE112393_PerCellAttributes.txt.gz','rt') as stream:
        header=[next(stream).rstrip() for _ in range(3)]
    report={'official_metadata_url':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112393/soft/GSE112393_family.soft.gz',
            'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'annotation_header':header,
            'batches':len(merged),'unresolved_aliases':unresolved,
            'note':'ST7/ST8 aliases unresolved; GEO has 1nDepleted1/2. GEO does not establish unique animal IDs for these batches. Do not call source-batch holdout leave-one-animal-out validation.'}
    (out/'atlas_metadata_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(merged[['source_batch','age','cell type enrichment','captured_cells']].to_string(index=False))


if __name__=='__main__':
    main()
