#!/usr/bin/env python3
"""Safely extract selected array members; never interpret arrays as counts."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def main():
    sheet=pd.read_csv(ROOT/'data/samples.tsv',sep='\t')
    manifest=[]
    for gse in ['GSE38096','GSE42601']:
        dest=ROOT/f'data/interim/arrays/{gse}'
        dest.mkdir(parents=True,exist_ok=True)
        selected=sheet[(sheet.gse==gse)&sheet.group.isin(['WT','KO'])].copy()
        archive=ROOT/f'data/raw/{gse}/{gse}_RAW.tar'
        names={}
        with tarfile.open(archive) as tar:
            for member in tar.getmembers():
                if not member.isfile():
                    continue
                if member.name=='GPL13912_old_annotations.txt.gz':
                    text=gzip.decompress(tar.extractfile(member).read()).decode().splitlines()
                    start=next(i for i,l in enumerate(text) if l.startswith('ID\t'))
                    table=pd.read_csv(io.StringIO('\n'.join(l for l in text[start:] if not l.startswith('!'))),sep='\t')
                    table.to_csv(dest/'probe_annotation.tsv',sep='\t',index=False)
                    continue
                match=re.search(r'GSM\d+',member.name)
                if not match or match.group() not in selected.sample_id.tolist():
                    continue
                name=Path(member.name).name
                content=gzip.decompress(tar.extractfile(member).read())
                target=dest/name.removesuffix('.gz')
                target.write_bytes(content)
                names[match.group()]=target.name
        if set(names)!=set(selected.sample_id):
            raise ValueError('Missing array samples')
        selected['raw_file']=selected.sample_id.map(names)
        selected['matrix_sample']=selected.sample_id
        selected.to_csv(dest/'samples.tsv',sep='\t',index=False)
        manifest.append({'gse':gse,'primary_samples':len(selected),'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
                         'excluded_HET':2 if gse=='GSE38096' else 0})
        print(gse,len(selected),'primary arrays',flush=True)
    (ROOT/'data/interim/arrays/manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    main()
