#!/usr/bin/env python3
"""Pinned mm10-era transcript mapping; compatibility must be checked by mass."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import urllib.request

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def main():
    dest=ROOT/'data/ref/annotation'
    dest.mkdir(parents=True,exist_ok=True)
    name='Mus_musculus.GRCm38.79.gtf.gz'
    path=dest/name
    url=f'https://ftp.ensembl.org/pub/release-79/gtf/mus_musculus/{name}'
    if not path.exists():
        with urllib.request.urlopen(url,timeout=50) as response:
            content=response.read()
        gzip.decompress(content)
        path.write_bytes(content)
    rows=[]
    with gzip.open(path,'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields=line.rstrip().split('\t')
            if len(fields)!=9 or fields[2]!='transcript':
                continue
            attrs=dict(re.findall(r'(\w+) "([^"]+)"',fields[8]))
            if 'gene_name' in attrs:
                rows.append((attrs['transcript_id'],attrs['gene_id'],attrs['gene_name']))
    mapping=pd.DataFrame(rows,columns=['transcript','ensembl_gene','symbol']).drop_duplicates()
    ambiguous=mapping.transcript.duplicated(keep=False)
    mapping[ambiguous].to_csv(dest/'ambiguous_transcript_mapping.tsv',sep='\t',index=False)
    mapping[~ambiguous].to_csv(dest/'ensembl79_transcript_mapping.tsv',sep='\t',index=False)
    notes={'url':url,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
           'unique_mapped_transcripts':int((~ambiguous).sum()),
           'note':'Original GEO specifies mm10 but not annotation release. Ensembl 79 is an explicit legacy compatibility reference, not asserted to be the original GTF. Check retained count mass and report this choice.'}
    (dest/'provenance.json').write_text(json.dumps(notes,indent=2)+'\n')
    print(notes)


if __name__=='__main__':
    main()
