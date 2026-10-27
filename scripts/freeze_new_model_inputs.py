#!/usr/bin/env python3
"""Reconcile primary-source retrieval evidence, then freeze new analysis inputs."""
import gzip
import json
from pathlib import Path
import pandas as pd
from discover_expansion import parse_soft
from freeze_replicability_inputs import digest
ROOT=Path(__file__).resolve().parents[1]
FIELDS={'title','source_name_ch1','organism_ch1','characteristics_ch1','library_strategy',
        'library_source','library_selection','data_processing'}


def main():
    out=ROOT/'results/replicability_discovery'
    discovery=json.loads((out/'discovery.json').read_text())
    inventory=json.loads((out/'metadata_inventory.json').read_text())
    assets=json.loads((out/'new_assets.json').read_text())['assets']
    paths=[]
    for row in discovery['request_provenance']:
        p=ROOT/row['cache_path']
        if p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Cached query source differs')
        paths.append(row['cache_path'])
    for row in inventory['sources']+assets:
        p=ROOT/row['path']
        if p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Downloaded source differs')
        paths.append(row['path'])
    for row in inventory['sources']:
        name=row['accession']+'_metadata.json';doc=json.loads((out/name).read_text())
        series,samples=parse_soft(gzip.decompress((ROOT/row['path']).read_bytes()).decode('utf-8'))
        simplified={gsm:{k.removeprefix('Sample_'):v for k,v in fields.items()
            if k.removeprefix('Sample_') in FIELDS or k.startswith('Sample_supplementary_file')} for gsm,fields in samples.items()}
        if doc['series']!=series or doc['samples']!=simplified or doc['sha256']!=row['sha256']:
            raise ValueError('Derived metadata does not reconcile to original SOFT')
        paths.append(str((out/name).relative_to(ROOT)))
    paths+=['results/replicability/candidate_stage_context.csv','results/replicability/primary_gene_summary.csv.gz',
        'results/replicability_discovery/discovery.json','results/replicability_discovery/metadata_inventory.json',
        'results/replicability_discovery/new_assets.json']
    rows=[{'path':p,'bytes':(ROOT/p).stat().st_size,'sha256':digest(ROOT/p)} for p in sorted(set(paths))]
    (ROOT/'data/new_model_inputs.json').write_text(json.dumps({'scope':'New source files, exact lookup responses, SOFT-reconciled metadata and fixed prior candidate tables. No frozen source guesses accepted.', 'assets':rows},indent=2)+'\n')
    review=json.loads((ROOT/'data/new_model_review.json').read_text())['candidates']
    pd.DataFrame(review).to_csv(out/'metadata_review.csv',index=False)
    compact=[]
    for r in discovery['candidates']:
        compact.append({k:r.get(k) for k in ['UID','accession','title','assay','samples','publication_date']}
                       |{'targeted_metadata_review':r['accession'] in {x['gse'] for x in review}})
    pd.DataFrame(compact).to_csv(out/'search_inventory.csv',index=False)
    audit={k:discovery[k] for k in ['retrieved_UTC','scope','endpoint','queries','deduplicated_UIDs','resolved_series',
        'HTTP_requests_this_run','public_unauthenticated','request_provenance','limits']}
    audit['metadata_records_reconciled_to_original_SOFT']=len(inventory['sources'])
    audit['metadata_HTTP_requests_this_run']=inventory['HTTP_requests_this_run']
    audit['metadata_sources']=inventory['sources']
    audit['identifier_resolution']='Entrez UID to returned ESummary accession field; no arithmetic conversion.'
    (out/'query_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('Frozen',len(rows),'inputs; reconciled',len(inventory['sources']),'original SOFT records')


if __name__=='__main__':main()
