#!/usr/bin/env python3
"""Retrieve only the frozen public assets through the launcher's local proxy."""
import json
import os
from pathlib import Path
import time
import urllib.request
from run_independent_contexts import ROOT,digest


def main():
    if not all(os.environ.get(k) for k in ['HTTP_PROXY','HTTPS_PROXY','http_proxy','https_proxy']):
        raise ValueError('Run through envs/python.sh with the local proxy configured; no direct fallback')
    assets=json.loads((ROOT/'data/independent_contexts_inputs.json').read_text())['assets']
    total=0
    for a in assets:
        p=ROOT/a['path']
        if p.exists():
            if p.stat().st_size!=a['bytes'] or digest(p)!=a['sha256']:raise ValueError('Cached asset differs: '+a['path'])
            continue
        if not a.get('url'):raise ValueError('Non-downloadable prerequisite missing: '+a['path'])
        if a['bytes']>100_000_000 or total+a['bytes']>250_000_000:raise ValueError('Download contract limit')
        p.parent.mkdir(parents=True,exist_ok=True);temporary=p.with_name(p.name+'.partial')
        try:
            with urllib.request.urlopen(a['url'],timeout=90) as response,temporary.open('wb') as f:
                size=0
                while block:=response.read(1024*1024):
                    size+=len(block)
                    if size>min(100_000_000,a['bytes']):raise ValueError('Unexpected asset size')
                    f.write(block)
            if size!=a['bytes'] or digest(temporary)!=a['sha256']:raise ValueError('Downloaded asset identity changed; investigate, do not refresh hashes')
            temporary.replace(p);total+=size;print(a['path'],size,flush=True);time.sleep(.6)
        finally:
            if temporary.exists():temporary.unlink()
    print('All frozen assets available and hash checked; newly downloaded bytes:',total)


if __name__=='__main__':main()
