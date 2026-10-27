#!/usr/bin/env python3
"""Bounded official GEO metadata retrieval for four prespecified candidates."""
import argparse
from datetime import datetime,timezone
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import time
import urllib.error
import urllib.request

from discover_expansion import parse_soft

ROOT=Path(__file__).resolve().parents[1]
CANDIDATES=['GSE155967','GSE107644','GSE109033','GSE83264']


def fetch(path,url,limit,download):
    if path.is_symlink():raise ValueError('Symlink cache refused')
    if not path.exists():
        if not download:raise FileNotFoundError('Missing cache; use --download')
        if not all(os.environ.get(k) for k in ['HTTP_PROXY','http_proxy','HTTPS_PROXY','https_proxy']):raise ValueError('Proxy required; no direct fallback')
        path.parent.mkdir(parents=True,exist_ok=True);partial=path.with_name(path.name+'.download-part')
        if partial.is_symlink():raise ValueError('Symlink partial refused')
        request=urllib.request.Request(url,headers={'User-Agent':'spermatogenic-arrest-plan-source-audit/1.0'})
        for attempt in range(2):
            try:
                with urllib.request.urlopen(request,timeout=45) as response,partial.open('wb') as out:
                    declared=response.headers.get('Content-Length');size=0;begun=time.monotonic()
                    if declared and int(declared)>limit:raise ValueError('Declared download exceeds cap')
                    for block in iter(lambda:response.read(1024**2),b''):
                        size+=len(block)
                        if size>limit or time.monotonic()-begun>180:raise ValueError('Download cap/deadline exceeded')
                        out.write(block)
                    if declared and size!=int(declared):raise ValueError('Incomplete download')
                partial.replace(path);break
            except urllib.error.HTTPError as exc:
                if exc.code not in [429,503] or attempt:raise
                time.sleep(2)
            except (urllib.error.URLError,TimeoutError) as exc:
                if attempt:raise
                print('Transient transfer failure; retry once through the same configured proxy:',type(exc).__name__,flush=True)
                time.sleep(2)
    if path.stat().st_size>limit:raise ValueError('Cached file exceeds cap')
    return {'path':str(path.relative_to(ROOT)),'url':url,'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--download',action='store_true');args=parser.parse_args()
    out=ROOT/'results/leptotene_followup';out.mkdir(parents=True,exist_ok=True)
    records=[];assets=[]
    for accession in CANDIDATES:
        url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/soft/{accession}_family.soft.gz'
        path=ROOT/f'data/raw/leptotene/metadata/{accession}_family.soft.gz'
        try:
            asset=fetch(path,url,16*1024**2,args.download)
            with gzip.open(path,'rb') as stream:data=stream.read(128*1024**2+1)
            if len(data)>128*1024**2:raise ValueError('Decompressed metadata exceeds cap')
            text=data.decode('utf-8');series,samples=parse_soft(text)
            if f'^SERIES = {accession}' not in text or len(series.get('Series_sample_id',[]))!=len(samples):raise ValueError('Series/sample count reconciliation failed')
            record={'accession':accession,'series':series,'samples':samples,'sample_count':len(samples),'metadata_asset':asset,
                'status':'metadata_retrieved_not_automatically_included'}
            assets.append(asset);records.append(record)
            (out/f'{accession}_metadata.json').write_text(json.dumps(record,indent=2)+'\n')
            print(accession,series.get('Series_title'),len(samples),'samples',asset['bytes'],'bytes',flush=True)
        except Exception as error:
            records.append({'accession':accession,'status':'failed_not_included','error':type(error).__name__+': '+str(error)})
            print(accession,'FAILED',type(error).__name__,str(error),flush=True)
        time.sleep(.4)
    report={'accessed_utc':datetime.now(timezone.utc).isoformat(),'scope':'Four prespecified candidate series; targeted, not exhaustive discovery.',
        'candidate_count':len(CANDIDATES),'metadata_retrieved':len(assets),'candidates':records,'assets':assets}
    (out/'discovery.json').write_text(json.dumps(report,indent=2)+'\n')
    if len(assets)!=len(CANDIDATES):raise RuntimeError('One or more candidate lookups failed; partial results retained explicitly')


if __name__=='__main__':main()
