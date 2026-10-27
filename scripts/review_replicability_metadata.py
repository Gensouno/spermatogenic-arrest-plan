#!/usr/bin/env python3
"""Fetch only the twelve contracted GEO family records, using configured proxy."""
import argparse
from datetime import datetime,timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import time
import urllib.request
from discover_expansion import parse_soft
ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--download',action='store_true');args=parser.parse_args()
    contract=json.loads((ROOT/'data/replicability_discovery_contract.json').read_text())
    accessions=contract['targeted_metadata_accessions']
    if len(accessions)>contract['limits']['maximum_new_SOFT_records'] or len(accessions)!=len(set(accessions)):
        raise ValueError('Metadata request budget differs')
    inventory=json.loads((ROOT/'results/replicability_discovery/discovery.json').read_text())
    if not set(accessions)<=set(r['accession'] for r in inventory['candidates']):raise ValueError('Metadata accession absent from search inventory')
    cache=ROOT/'data/raw/replicability_discovery';cache.mkdir(parents=True,exist_ok=True)
    out=ROOT/'results/replicability_discovery';out.mkdir(parents=True,exist_ok=True)
    records=[];requested=0
    for acc in accessions:
        if not re.fullmatch(r'GSE[0-9]+',acc):raise ValueError('Invalid accession')
        url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz'
        path=cache/(acc+'_family.soft.gz')
        if not path.exists():
            if not args.download:raise FileNotFoundError('Explicit --download required for uncached SOFT')
            for var in ['HTTP_PROXY','HTTPS_PROXY','http_proxy','https_proxy']:
                if not os.environ.get(var):raise ValueError('Configured proxy required before retrieval')
            time.sleep(.55)
            with urllib.request.urlopen(url,timeout=45) as response:data=response.read()
            text=gzip.decompress(data).decode('utf-8')
            if '^SERIES = '+acc not in text:raise ValueError('Unexpected GEO response')
            path.write_bytes(data);requested+=1
        data=path.read_bytes();series,samples=parse_soft(gzip.decompress(data).decode('utf-8'))
        simplified={gsm:{k.removeprefix('Sample_'):v for k,v in fields.items()
            if k in ['Sample_title','Sample_source_name_ch1','Sample_organism_ch1','Sample_characteristics_ch1',
                     'Sample_library_strategy','Sample_library_source','Sample_library_selection','Sample_data_processing']
                or k.startswith('Sample_supplementary_file')} for gsm,fields in samples.items()}
        record={'accession':acc,'source_url':url,'path':str(path.relative_to(ROOT)),
            'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'series':series,'samples':simplified,
            'status':'Sample-level source evidence; eligibility not automatically established'}
        (out/(acc+'_metadata.json')).write_text(json.dumps(record,indent=2)+'\n');records.append(record)
        print(acc,len(samples),series.get('Series_supplementary_file',[]),flush=True)
    (out/'metadata_inventory.json').write_text(json.dumps({'retrieved_UTC':datetime.now(timezone.utc).isoformat(),
        'contract_accessions':accessions,'resolved_records':len(records),'HTTP_requests_this_run':requested,
        'sources':[{k:r[k] for k in ['accession','source_url','path','bytes','sha256']} for r in records]},indent=2)+'\n')


if __name__=='__main__':main()
