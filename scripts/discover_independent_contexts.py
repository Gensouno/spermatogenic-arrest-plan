from pathlib import Path
import hashlib,json,os,re,time,urllib.parse,urllib.request
r=Path.cwd();c=json.loads((r/'data/independent_discovery_contract.json').read_text());out=r/'results/independent_discovery';raw=r/'data/raw/independent_discovery';records=[];calls=0

def request(endpoint,params):
 global calls
 url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+endpoint+'.fcgi?'+urllib.parse.urlencode(params)
 p=raw/(hashlib.sha256(url.encode()).hexdigest()+'.json')
 if not p.exists():
  assert all(os.environ.get(k) for k in ['HTTP_PROXY','HTTPS_PROXY','http_proxy','https_proxy'])
  calls+=1;assert calls<=c['limits']['maximum_API_requests'];time.sleep(.55)
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'public-research-source-audit/1.0'}),timeout=45) as response:data=response.read(4*1024*1024)
  d=json.loads(data);assert 'error' not in d;p.write_bytes(data)
 data=p.read_bytes();d=json.loads(data);records.append(dict(url=url,path=str(p.relative_to(r)),bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
 return d
queries=[];allids=set()
for query in c['queries']:
 count=request('esearch',dict(db='gds',term=query,retmode='json',retmax=0))['esearchresult'];n=int(count['count']);print('Count:',n,query,flush=True)
 assert n<=c['limits']['maximum_series_per_query'],'Query exceeds contract bound; partition explicitly'
 d=request('esearch',dict(db='gds',term=query,retmode='json',retmax=500,retstart=0,sort='Accession'))['esearchresult']
 ids=d['idlist'];assert len(ids)==len(set(ids))==n and all(re.fullmatch(r'\d+',v) for v in ids);allids.update(ids)
 queries.append(dict(query=query,server_count=n,retrieved=len(ids),retstart=0,retmax=500,query_translation=d.get('querytranslation')))
rows=[];ids=sorted(allids,key=int)
for start in range(0,len(ids),40):
 batch=ids[start:start+40];d=request('esummary',dict(db='gds',id=','.join(batch),retmode='json'))['result'];assert set(d['uids'])==set(batch)
 for uid in batch:
  v=d[uid];acc=v['accession'];assert re.fullmatch(r'GSE\d+',acc)
  rows.append(dict(UID=uid,accession=acc,title=v.get('title'),summary=v.get('summary'),assay=v.get('gdstype'),samples=v.get('n_samples'),publication_date=v.get('pdat'),pubmedids=v.get('pubmedids')))
(out/'discovery.json').write_text(json.dumps(dict(accessed_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),scope=c['scope'],queries=queries,unique_UIDs=len(ids),resolved_series=len(rows),request_provenance=records,HTTP_requests_this_run=calls,candidates=rows),indent=2)+'\n')
import pandas as pd
pd.DataFrame(rows).drop(columns=['summary','pubmedids']).to_csv(out/'search_inventory.csv',index=False)
print('Resolved unique series:',len(rows),'HTTP calls:',calls,flush=True)
