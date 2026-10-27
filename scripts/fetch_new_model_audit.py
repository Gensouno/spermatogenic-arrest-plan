#!/usr/bin/env python3
"""Bounded proxied acquisition of two explicit processed-data assets."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import urllib.request
ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--download',action='store_true');args=parser.parse_args()
    spec=json.loads((ROOT/'data/new_model_audit_contract.json').read_text());rows=[]
    for row in spec['assets']:
        path=ROOT/row['path'];path.parent.mkdir(parents=True,exist_ok=True)
        if not path.exists():
            if not args.download:raise FileNotFoundError('Explicit --download required for uncached asset')
            for var in ['HTTP_PROXY','HTTPS_PROXY','http_proxy','https_proxy']:
                if not os.environ.get(var):raise ValueError('Configured proxy required before retrieval')
            chunks=[];size=0
            with urllib.request.urlopen(row['url'],timeout=60) as response:
                while chunk:=response.read(1024*1024):
                    size+=len(chunk)
                    if size>row['maximum_bytes']:raise ValueError('Asset exceeds bounded retrieval size')
                    chunks.append(chunk)
            path.write_bytes(b''.join(chunks))
        size=path.stat().st_size
        if not 0<size<=row['maximum_bytes']:raise ValueError('Invalid cached asset size')
        rows.append({**row,'bytes':size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
        print(row['gse'],size,rows[-1]['sha256'],flush=True)
    out=ROOT/'results/replicability_discovery';out.mkdir(parents=True,exist_ok=True)
    (out/'new_assets.json').write_text(json.dumps({'retrieved_UTC':datetime.now(timezone.utc).isoformat(),
        'assets':rows,'URL_conversion':'GEO ftp supplementary links converted to same-host HTTPS; no data-content conversion.'},indent=2)+'\n')


if __name__=='__main__':main()
