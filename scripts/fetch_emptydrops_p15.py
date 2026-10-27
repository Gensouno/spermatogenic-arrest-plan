#!/usr/bin/env python3
"""Retrieve one complete P15 EmptyDrops library from its official matrix."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from emptydrops_range_client import MatrixRanges,digest
from fetch_labeled_reference import file_entries

ROOT=Path(__file__).resolve().parents[1]
BASE='https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/946/E-MTAB-6946/Files/'


def verify_frozen(row,download=False):
    """Restore a missing file only when explicitly requested; never replace mismatch."""
    from check_project import safe_path
    path=safe_path(ROOT,row['path'])
    if not path.exists():
        if not download:raise FileNotFoundError('Frozen input missing; use --download: '+row['path'])
        import os,urllib.request
        if not all(os.environ.get(k) for k in ['HTTP_PROXY','http_proxy','HTTPS_PROXY','https_proxy']):raise ValueError('Proxy required')
        path.parent.mkdir(parents=True,exist_ok=True)
        partial=path.with_name(path.name+'.restore-part')
        if partial.is_symlink():raise ValueError('Symlink partial refused')
        if 'range' in row:
            remote=MatrixRanges(row['url'],2224359116,etag=row['etag'])
            w=row['range'];remote.download(partial,w['start_byte'],w['end_byte_exclusive']-1)
        else:
            with urllib.request.urlopen(row['url'],timeout=40) as response,partial.open('wb') as stream:
                size=0
                for block in iter(lambda:response.read(1024**2),b''):
                    size+=len(block)
                    if size>row['bytes']:raise ValueError('Restored source oversized')
                    stream.write(block)
        if partial.stat().st_size!=row['bytes'] or digest(partial)!=row['sha256']:raise ValueError('Restored source differs; partial retained')
        partial.replace(path)
    if path.stat().st_size!=row['bytes'] or digest(path)!=row['sha256']:raise ValueError('Frozen EmptyDrops input differs: '+row['path'])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--download',action='store_true');args=parser.parse_args()
    out=ROOT/'data/raw/emptydrops_p15';out.mkdir(parents=True,exist_ok=True)
    spec=ROOT/'data/emptydrops_p15_assets.json'
    if spec.exists():
        frozen=json.loads(spec.read_text())
        for row in frozen['assets']:
            verify_frozen(row,args.download)
        print('Verified',len(frozen['assets']),'EmptyDrops input files; existing valid files reused');return
    if not args.download:raise FileNotFoundError('EmptyDrops P15 cache not frozen; use --download')
    # Reuse previously hash-frozen official study/annotation files, not new labels.
    prior=json.loads((ROOT/'data/labeled_reference_assets.json').read_text())['assets']
    prior+=json.loads((ROOT/'data/emptydrops_label_assets.json').read_text())['assets']
    needed={'study.json','cell_metadata.txt','cell_metadata_emptyDrops.txt','E-MTAB-6946.sdrf.txt'}
    records=[]
    for row in prior:
        if Path(row['path']).name in needed:
            p=ROOT/row['path']
            if digest(p)!=row['sha256'] or p.stat().st_size!=row['bytes']:raise ValueError('Prior source input changed')
            records.append(row)
    if len(records)!=4:raise ValueError('Required metadata inventory differs')
    declared={r['path']:r for r in file_entries(json.loads((ROOT/'data/raw/labeled_reference/study.json').read_text()))}
    import urllib.request,os
    if not (os.environ.get('HTTPS_PROXY') and os.environ.get('https_proxy')):raise ValueError('Proxy required')
    name='genes_emptyDrops.tsv';path=out/name
    if not path.exists():
        with urllib.request.urlopen(BASE+name,timeout=40) as response:data=response.read(2*1024**2)
        if len(data)!=declared[name]['size']:raise ValueError('Gene table declared size differs')
        path.write_bytes(data)
    if path.is_symlink() or path.stat().st_size!=declared[name]['size']:raise ValueError('Unexpected gene table')
    records.append({'path':str(path.relative_to(ROOT)),'url':BASE+name,'bytes':path.stat().st_size,'sha256':digest(path)})
    cells=pd.read_csv(ROOT/'data/raw/labeled_reference/cell_metadata_emptyDrops.txt',sep=r'\s+',dtype={'P15Clusters':'string'},low_memory=False)
    if cells.index.tolist()!=list(range(1,len(cells)+1)):raise ValueError('Unexpected author row order')
    selected=cells.index[(cells.Library=='do18195')&(cells.Sample=='P15')].to_numpy()
    if len(selected)!=13168 or not np.array_equal(selected,np.arange(selected[0],selected[-1]+1)):raise ValueError('P15 library not complete/contiguous')
    cumulative=np.r_[0,cells.total_features_by_counts.to_numpy(dtype=np.int64).cumsum()]
    genes=pd.read_csv(path,sep='\t')
    remote=MatrixRanges(BASE+'raw_counts_emptyDrops.mtx',declared['raw_counts_emptyDrops.mtx']['size'])
    with remote.request(0,4095) as response:head=response.read()
    header=head.split(b'\n',2)[:2]
    if header[0]!=b'%%MatrixMarket matrix coordinate integer general':raise ValueError('Unexpected MatrixMarket format')
    rows,cols,nnz=map(int,header[1].split())
    if (rows,cols,nnz)!=(len(genes),len(cells),int(cumulative[-1])):raise ValueError('Header/author metadata dimensions differ')
    header_bytes=sum(len(x)+1 for x in header)
    start=remote.boundary(int(selected[0]),cumulative,header_bytes)
    end=remote.boundary(int(selected[-1])+1,cumulative,header_bytes)
    if end-start>500*1024**2:raise ValueError('Window exceeds prespecified 500 MiB cap')
    window={'library':'do18195','first_column':int(selected[0]),'last_column':int(selected[-1]),'start_byte':start,
        'end_byte_exclusive':end,'cells':len(selected),'nonzeros':int(cumulative[selected[-1]]-cumulative[selected[0]-1])}
    plan={'source_matrix':{'url':remote.url,'bytes':remote.size,'etag':remote.etag,'header':[x.decode() for x in header],
        'full_matrix_sha256_available':False},'window':window}
    logs=ROOT/'logs/emptydrops_p15';logs.mkdir(parents=True,exist_ok=True)
    (logs/'range_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    print('Resolved complete P15 window:',json.dumps(window),flush=True)
    path=out/'do18195_emptyDrops.triplets'
    result=remote.download(path,start,end-1)
    records.append({'path':str(path.relative_to(ROOT)),'url':remote.url,'range':window,'etag':remote.etag,
        'bytes':result['bytes'],'sha256':result['sha256']})
    report={'schema':1,'retrieved_utc':datetime.now(timezone.utc).isoformat(),
        'scope':'One complete P15 author-called EmptyDrops library; not the complete atlas or raw sequencing reads.',
        **plan,'range_requests_this_run':len(remote.calls),'probe_and_transfer_ranges':remote.calls,
        'verified_transfer_chunks':result['chunks'],'assets':records}
    spec.write_text(json.dumps(report,indent=2)+'\n')
    print('Frozen',len(records),'inputs; P15 triplet bytes:',result['bytes'],flush=True)


if __name__=='__main__':main()
