#!/usr/bin/env python3
"""Bounded version-checked HTTP ranges for one official sparse matrix."""
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.request

import numpy as np


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for b in iter(lambda:stream.read(1024**2),b''):h.update(b)
    return h.hexdigest()


class MatrixRanges:
    def __init__(self,url,size,etag=None):
        if not url.startswith('https://ftp.ebi.ac.uk/biostudies/') or size<=0:raise ValueError('Unexpected matrix source')
        self.url=url;self.size=size;self.etag=etag;self.calls=[]

    def request(self,start,end):
        if not all(os.environ.get(k) for k in ['HTTP_PROXY','http_proxy','HTTPS_PROXY','https_proxy']):raise ValueError('Explicit proxy required')
        if not 0<=start<=end<self.size or len(self.calls)>=85:raise ValueError('Invalid range/request budget')
        headers={'Range':f'bytes={start}-{end}','User-Agent':'spermatogenic-arrest-plan-source-audit/1.0'}
        if self.etag:headers['If-Match']=self.etag
        response=urllib.request.urlopen(urllib.request.Request(self.url,headers=headers),timeout=45)
        tag=response.headers.get('ETag')
        if (response.status!=206 or response.headers.get('Content-Range')!=f'bytes {start}-{end}/{self.size}'
                or not tag or tag.startswith('W/') or (self.etag and self.etag!=tag)):
            response.close();raise ValueError('Source version or exact range response differs')
        self.etag=tag;self.calls.append({'start':start,'end':end,'etag':tag})
        return response

    def boundary(self,column,cumulative,header_bytes):
        if column==1:return header_bytes
        if column==len(cumulative):return self.size
        if not 1<column<len(cumulative):raise ValueError('Column outside metadata')
        lo=header_bytes;hi=self.size
        guess=int(cumulative[column-1]/cumulative[-1]*(self.size-header_bytes))+header_bytes
        for _ in range(20):
            start=max(lo,min(guess,hi-1));end=min(start+131071,self.size-1)
            with self.request(start,end) as r:block=r.read()
            if len(block)!=end-start+1:raise ValueError('Truncated boundary probe')
            skip=block.index(b'\n')+1;body=block[skip:block.rfind(b'\n')+1]
            lines=body.splitlines(keepends=True)
            a=np.fromstring(body.decode('ascii'),sep=' ',dtype=np.int64).reshape(-1,3)
            columns=a[:,1]
            if len(a)!=len(lines) or (np.diff(columns)<0).any():raise ValueError('Unexpected matrix ordering')
            crossing=np.flatnonzero((columns[:-1]<column)&(columns[1:]>=column))
            if len(crossing):
                i=int(crossing[0])+1
                if columns[i]!=column or columns[i-1]!=column-1:raise ValueError('Missing boundary column')
                return start+skip+sum(map(len,lines[:i]))
            first=int(columns[0]);last=int(columns[-1])
            if first>=column:hi=start+skip
            elif last<column:lo=start+skip+len(body)
            else:raise ValueError('Ambiguous boundary')
            delta=int((cumulative[column-1]-cumulative[first-1])*(self.size-header_bytes)/cumulative[-1])
            candidate=start+skip+delta-65536
            if not lo<candidate<hi or candidate==guess:candidate=(lo+hi)//2
            guess=candidate
        raise ValueError('Exact column boundary not resolved')

    def download(self,path,start,end):
        """Resume verified complete 8 MiB chunks; assemble no partial into final."""
        if not 0<=start<=end<self.size:raise ValueError('Invalid download range')
        path=Path(path);chunks=path.with_name(path.name+'.parts')
        if chunks.is_symlink():raise ValueError('Symlink chunk directory refused')
        chunks.mkdir(parents=True,exist_ok=True)
        index_path=chunks/'index.json'
        if index_path.is_symlink():raise ValueError('Symlink chunk index refused')
        index={'url':self.url,'matrix_bytes':self.size,'etag':self.etag,'start':start,'end':end,'chunks':{}}
        if index_path.exists():
            old=json.loads(index_path.read_text())
            if any(old[k]!=index[k] for k in ['url','matrix_bytes','etag','start','end']):raise ValueError('Existing chunk source differs')
            index=old
        begun=time.monotonic();total=0;paths=[]
        for number,offset in enumerate(range(start,end+1,8*1024**2)):
            stop=min(offset+8*1024**2-1,end);p=chunks/f'{number:04d}.bin';key=p.name
            if p.is_symlink():raise ValueError('Symlink chunk refused')
            if p.exists():
                row=index['chunks'].get(key)
                if not row or p.stat().st_size!=stop-offset+1 or digest(p)!=row['sha256']:raise ValueError('Unverified/different cached chunk')
            else:
                partial=p.with_suffix('.download-part')
                if partial.is_symlink():raise ValueError('Symlink partial refused')
                with self.request(offset,stop) as response,partial.open('wb') as out:
                    size=0
                    for block in iter(lambda:response.read(1024**2),b''):
                        size+=len(block)
                        if size>stop-offset+1:raise ValueError('Oversized chunk')
                        out.write(block)
                if size!=stop-offset+1:raise ValueError('Truncated chunk; partial retained')
                partial.replace(p)
                index['chunks'][key]={'start':offset,'end':stop,'bytes':size,'sha256':digest(p)}
                index_path.write_text(json.dumps(index,indent=2)+'\n')
            total+=p.stat().st_size;paths.append(p)
            print('Verified matrix window',round(total/1024**2,1),'MiB;',round(time.monotonic()-begun,1),'seconds',flush=True)
            if time.monotonic()-begun>1200:raise TimeoutError('Twenty-minute transfer bound; verified chunks retained')
        partial=path.with_name(path.name+'.download-part')
        if path.is_symlink() or partial.is_symlink():raise ValueError('Symlink final/partial refused')
        with partial.open('wb') as stream:
            for p in paths:
                with p.open('rb') as source:
                    for block in iter(lambda:source.read(1024**2),b''):stream.write(block)
        if partial.stat().st_size!=end-start+1:raise ValueError('Assembled size differs')
        partial.replace(path)
        return {'bytes':total,'sha256':digest(path),'chunks':len(paths)}
