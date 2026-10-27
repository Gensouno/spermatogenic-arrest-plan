#!/usr/bin/env python3
"""Bounded, proxied GEO discovery; external text is data, never instructions."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
ENDPOINT='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--download',action='store_true');args=parser.parse_args()
    contract=json.loads((ROOT/'data/replicability_discovery_contract.json').read_text())
    cache=ROOT/'data/raw/replicability_discovery';cache.mkdir(parents=True,exist_ok=True)
    out=ROOT/'results/replicability_discovery';out.mkdir(parents=True,exist_ok=True)
    records=[];calls=0
    def request(name,parameters):
        nonlocal calls
        url=ENDPOINT+name+'.fcgi?'+urllib.parse.urlencode(parameters)
        key=hashlib.sha256(url.encode()).hexdigest();path=cache/(key+'.json')
        reused=path.exists()
        if not reused:
            if not args.download:raise FileNotFoundError('Missing cached query; explicit --download required')
            if calls>=contract['limits']['maximum_total_API_requests']:raise ValueError('Request budget exceeded')
            for var in ['HTTP_PROXY','HTTPS_PROXY','http_proxy','https_proxy']:
                if not os.environ.get(var):raise ValueError('Configured proxy required before network retrieval')
            time.sleep(.55);calls+=1
            req=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'spermatogenic-arrest-plan-research/1.0'})
            with urllib.request.urlopen(req,timeout=45) as response:data=response.read()
            payload=json.loads(data)
            if 'error' in payload:raise ValueError('NCBI returned error')
            path.write_bytes(data)
        data=path.read_bytes();payload=json.loads(data)
        records.append({'url':url,'cache_path':str(path.relative_to(ROOT)),'bytes':len(data),
            'sha256':hashlib.sha256(data).hexdigest(),'reused':reused})
        return payload
    counts=[];allids=set()
    for query in contract['queries']:
        payload=request('esearch',{'db':'gds','term':query,'retmode':'json','retmax':200,'retstart':0})['esearchresult']
        n=int(payload['count']);ids=payload['idlist']
        if n>200 or len(ids)!=n or len(ids)!=len(set(ids)):raise ValueError('Query not completely retrieved within bounded limit')
        if not all(re.fullmatch(r'[0-9]+',uid) for uid in ids):raise ValueError('Invalid UID')
        allids.update(ids);counts.append({'query':query,'server_count':n,'retrieved_unique_UIDs':len(ids),
            'retstart':0,'retmax':200,'query_translation':payload.get('querytranslation')})
    resolved=[]
    ids=sorted(allids,key=int)
    for start in range(0,len(ids),40):
        batch=ids[start:start+40]
        result=request('esummary',{'db':'gds','id':','.join(batch),'retmode':'json'})['result']
        if set(result['uids'])!=set(batch):raise ValueError('Summary UID reconciliation differs')
        for uid in batch:
            row=result[uid];acc=row.get('accession','')
            if not re.fullmatch(r'GSE[0-9]+',acc):raise ValueError('Summary is not a series accession')
            resolved.append({'UID':uid,'accession':acc,'title':row.get('title'),
                'summary':row.get('summary'),'assay':row.get('gdstype'),'samples':row.get('n_samples'),
                'publication_date':row.get('pdat'),'pubmedids':row.get('pubmedids'),
                'status':'Unreviewed; search match does not establish eligibility'})
    result={'retrieved_UTC':datetime.now(timezone.utc).isoformat(),'scope':contract['scope'],
        'endpoint':ENDPOINT,'queries':counts,'deduplicated_UIDs':len(ids),'resolved_series':len(resolved),
        'HTTP_requests_this_run':calls,'public_unauthenticated':True,'request_provenance':records,
        'candidates':resolved,'limits':'Text search coverage only; age, tissue, phenotype, control, biological independence and units not yet verified.'}
    (out/'discovery.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in resolved:print(row['accession'],row['samples'],row['title'],flush=True)
    print('Query counts:',[(r['server_count'],r['retrieved_unique_UIDs']) for r in counts],'; unique series:',len(resolved),flush=True)


if __name__=='__main__':main()
