#!/usr/bin/env python3
"""Read three prespecified complete library windows from the official sorted matrix."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import re
import time
import urllib.request

import numpy as np
import pandas as pd
from fetch_labeled_reference import CACHE,ROOT,AUTHOR_COMMIT,digest,fetch

URL='https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/946/E-MTAB-6946/Files/raw_counts.mtx'
SIZE=3012455192
LIBRARIES=['do17815','do26386','do18195']


class RemoteMatrix:
    def __init__(self):self.etag=None;self.calls=[]

    def request(self,start,end):
        if len(self.calls)>=95:raise ValueError('Range request budget exhausted')
        headers={'Range':f'bytes={start}-{end}','User-Agent':'spermatogenic-arrest-plan-source-audit/1.0'}
        if self.etag:headers['If-Match']=self.etag
        r=urllib.request.urlopen(urllib.request.Request(URL,headers=headers),timeout=40)
        expected=f'bytes {start}-{end}/{SIZE}'
        if r.status!=206 or r.headers.get('Content-Range')!=expected:
            r.close();raise ValueError('Server did not honor exact range')
        tag=r.headers.get('ETag')
        if not tag or tag.startswith('W/') or (self.etag and self.etag!=tag):r.close();raise ValueError('Source ETag missing/weak/changed')
        self.etag=tag;self.calls.append({'start':start,'end':end,'etag':tag})
        return r

    def boundary(self,column,cumulative,header_bytes):
        if column==1:return header_bytes
        lo=header_bytes;hi=SIZE
        guess=int(cumulative[column-1]/cumulative[-1]*(SIZE-header_bytes))+header_bytes
        for attempt in range(20):
            start=max(lo,min(guess,hi-1));end=min(start+131071,SIZE-1)
            with self.request(start,end) as r:block=r.read()
            if len(block)!=end-start+1:raise ValueError('Incomplete boundary probe')
            skip=block.index(b'\n')+1
            body=block[skip:block.rfind(b'\n')+1]
            lines=body.splitlines(keepends=True)
            a=np.fromstring(body.decode(),sep=' ',dtype=np.int64).reshape(-1,3)
            cols=a[:,1]
            crossing=np.flatnonzero((cols[:-1]<column)&(cols[1:]>=column))
            if len(crossing):
                i=int(crossing[0])+1
                if cols[i]!=column or cols[i-1]!=column-1:raise ValueError('Missing matrix column at boundary')
                return start+skip+sum(map(len,lines[:i]))
            first=int(cols[0]);last=int(cols[-1])
            if first>=column:hi=start+skip
            elif last<column:lo=start+skip+len(body)
            else:raise ValueError('Ambiguous column boundary')
            # Use known author nnz-per-cell to accelerate an otherwise bounded
            # byte search; exact transition checks above, not this guess, decide.
            delta=int((cumulative[column-1]-cumulative[first-1])*(SIZE-header_bytes)/cumulative[-1])
            candidate=start+skip+delta-65536
            if not lo<candidate<hi or candidate==guess:candidate=(lo+hi)//2
            guess=candidate
        raise ValueError('Could not resolve exact column boundary')


def restore_window(row,remote):
    if not (os.environ.get('HTTPS_PROXY') and os.environ.get('https_proxy')):raise ValueError('Proxy required')
    window=row['range'];library=window['library']
    if library not in LIBRARIES or row['url']!=URL:raise ValueError('Unexpected frozen window source')
    path=CACHE/(library+'.triplets');partial=path.with_name(path.name+'.download-part')
    if path.is_symlink() or partial.is_symlink():raise ValueError('Symlink refused')
    remote.etag=row['etag'];size=0;started=time.monotonic()
    with remote.request(window['start_byte'],window['end_byte_exclusive']-1) as response,partial.open('wb') as out:
        for block in iter(lambda:response.read(1024**2),b''):
            size+=len(block)
            if size>row['bytes'] or time.monotonic()-started>480:raise ValueError('Window exceeds transfer bound')
            out.write(block)
    if size!=row['bytes'] or digest(partial)!=row['sha256']:raise ValueError('Frozen window content differs; partial retained')
    partial.replace(path)
    print('Restored frozen window',library,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--download',action='store_true');args=p.parse_args()
    frozen_path=ROOT/'data/labeled_reference_assets.json'
    if frozen_path.exists():
        frozen=json.loads(frozen_path.read_text())
        remote=RemoteMatrix()
        for r in frozen['assets']:
            path=ROOT/r['path']
            if not path.exists():
                if not args.download:raise FileNotFoundError('Missing frozen cache; use --download: '+r['path'])
                CACHE.mkdir(parents=True,exist_ok=True)
                if 'range' in r:restore_window(r,remote)
                else:fetch(path.name,r['url'],15*1024**2,True)
            if path.is_symlink() or path.stat().st_size!=r['bytes'] or digest(path)!=r['sha256']:raise ValueError('Frozen source differs')
        print('Verified',len(frozen['assets']),'frozen files');return
    if not args.download:raise FileNotFoundError('No frozen range cache; use --download')
    if not (os.environ.get('HTTPS_PROXY') and os.environ.get('https_proxy')):raise ValueError('Proxy required')
    files=[('study.json','https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-6946'),
        ('author_README.md',f'https://raw.githubusercontent.com/MarioniLab/Spermatogenesis2018/{AUTHOR_COMMIT}/README.md')]
    files += [(n,'https://www.ebi.ac.uk/biostudies/files/E-MTAB-6946/'+n) for n in ['cell_metadata.txt','genes.tsv','E-MTAB-6946.sdrf.txt']]
    assets=[fetch(n,u,15*1024**2,True) for n,u in files]
    cells=pd.read_csv(CACHE/'cell_metadata.txt',sep=r'\s+')
    cumulative=np.r_[0,np.cumsum(cells.total_features_by_counts.to_numpy(dtype=np.int64))]
    remote=RemoteMatrix()
    with remote.request(0,255) as r:first=r.read()
    header=first.split(b'\n',2)[:2]
    if header!=[b'%%MatrixMarket matrix coordinate integer general',b'33226 53510 220823435']:
        raise ValueError('Unexpected matrix header')
    header_bytes=sum(len(x)+1 for x in header)
    if cumulative[-1]!=220823435:raise ValueError('Author total detected features differ from matrix nnz')
    windows=[]
    for library in LIBRARIES:
        selected=cells.index[cells.Library==library].to_numpy()
        if not np.array_equal(selected,np.arange(selected[0],selected[-1]+1)):raise ValueError('Library not contiguous')
        start=remote.boundary(int(selected[0]),cumulative,header_bytes)
        end=remote.boundary(int(selected[-1])+1,cumulative,header_bytes)
        windows.append({'library':library,'first_column':int(selected[0]),'last_column':int(selected[-1]),
            'start_byte':start,'end_byte_exclusive':end,'cells':len(selected),
            'nonzeros':int(cumulative[selected[-1]]-cumulative[selected[0]-1])})
        print('Resolved',library,windows[-1],flush=True)
    for window in windows:
        name=window['library']+'.triplets';path=CACHE/name;partial=path.with_name(name+'.download-part')
        if path.is_symlink() or partial.is_symlink():raise ValueError('Symlink refused')
        size=0;started=time.monotonic()
        with remote.request(window['start_byte'],window['end_byte_exclusive']-1) as r,partial.open('wb') as out:
            for block in iter(lambda:r.read(1024**2),b''):
                out.write(block);size+=len(block)
                if time.monotonic()-started>480:raise TimeoutError('Eight-minute window transfer limit')
                if size%(25*1024**2)==0:print(name,size//1024**2,'MiB',round(time.monotonic()-started,1),'seconds',flush=True)
        if size!=window['end_byte_exclusive']-window['start_byte']:raise ValueError('Incomplete library transfer')
        partial.replace(path)
        assets.append({'path':str(path.relative_to(ROOT)),'url':URL,'bytes':size,'sha256':digest(path),
            'range':window,'etag':remote.etag})
        print('Fetched complete library',window['library'],flush=True)
    report={'schema':1,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'author_commit':AUTHOR_COMMIT,
        'scope':'Three complete prespecified normal libraries from original column-sorted matrix; not the complete atlas or all raw sequencing.',
        'source_matrix':{'url':URL,'bytes':SIZE,'etag':remote.etag,'header':[s.decode() for s in header],
            'full_matrix_sha256_available':False},'author_code_executed':False,'range_requests':len(remote.calls),
        'probe_and_transfer_ranges':remote.calls,'assets':assets}
    frozen_path.write_text(json.dumps(report,indent=2)+'\n')
    print('Frozen',len(assets),'files;',sum(r['bytes'] for r in assets),'bytes;',len(remote.calls),'range requests',flush=True)


if __name__=='__main__':main()
