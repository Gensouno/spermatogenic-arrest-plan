#!/usr/bin/env python3
"""Retrieve the frozen successful source payloads via the configured proxy only."""
import json
import os
import time
import urllib.request
from run_stage_specificity import ROOT,digest


def main():
    if not all(os.environ.get(k) for k in ['HTTP_PROXY','HTTPS_PROXY','http_proxy','https_proxy']):raise ValueError('Proxy configuration required; no direct fallback')
    assets=json.loads((ROOT/'data/stage_specificity_inputs.json').read_text())['assets']
    for asset in assets:
        if asset.get('role')!='source_refresh':continue
        path=ROOT/asset['path']
        if not path.exists():
            with urllib.request.urlopen(asset['url'],timeout=45) as response:data=response.read(4_000_001)
            if len(data)>4_000_000 or len(data)!=asset['bytes']:raise ValueError('Source size differs')
            import hashlib
            if hashlib.sha256(data).hexdigest()!=asset['sha256']:raise ValueError('Source changed; investigate without refreshing the frozen identity')
            path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data);time.sleep(.6)
        if path.stat().st_size!=asset['bytes'] or digest(path)!=asset['sha256']:raise ValueError('Cached source identity differs')
    print('Frozen source payloads available; historic failed lookups remain explicit in data/stage_source_retrieval.json')


if __name__=='__main__':main()
