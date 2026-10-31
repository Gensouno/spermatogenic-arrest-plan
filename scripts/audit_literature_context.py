#!/usr/bin/env python3
"""Audit the dated private retrieval inventory without querying live databases."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from run_independent_contexts import digest

ROOT=Path(__file__).resolve().parents[1]


def audit(root):
    summary=json.loads((root/'research/literature-context-summary-2026-10-31.json').read_text())
    raw=root/'data/raw/literature_context_20261031';union=set();payloads=failures=0
    for key,d in summary['queries'].items():
        source=json.loads((raw/(key+'_ids.json')).read_text())['esearchresult']
        count=json.loads((raw/(key+'_count.json')).read_text())['esearchresult']['count']
        ids=d['ids']
        if len(ids)!=len(set(ids)) or len(ids)!=d['retrieved'] or len(ids)!=int(count):raise ValueError('Query count/ID reconciliation')
        if ids!=source['idlist'] or int(source['count'])!=d['server_count']:raise ValueError('Query source identity')
        union.update(ids)
    seen=[]
    for p in raw.glob('abstracts_*.xml'):
        seen.extend(e.text for e in ET.parse(p).findall('./PubmedArticle/MedlineCitation/PMID'))
    if len(seen)!=len(set(seen)) or set(seen)!=union or len(seen)!=summary['retrieved_abstract_records']:raise ValueError('Abstract identity reconciliation')
    records=pd.read_csv(raw/'records.tsv',sep='\t')
    if set(records.PMID.astype(str))!=union or len(records)!=len(union):raise ValueError('Parsed record identity')
    for name in ['literature-context-retrieval-2026-10-31.json','literature-followup-retrieval-2026-10-31.json','paired-cytology-candidate-retrieval.json']:
        for r in json.loads((root/'research'/name).read_text()):
            if r['status']!='success':failures+=1;continue
            p=root/r['path']
            if digest(p)!=r['sha256'] or p.stat().st_size!=r['bytes']:raise ValueError('Payload identity changed')
            payloads+=1
    candidate=json.loads((root/'data/paired_cytology_candidate.json').read_text())
    for spec in candidate['QC']:
        d=pd.read_csv(raw/spec['file'],sep='\t');n=d.select_dtypes('number');x=n.to_numpy()
        if len(d)!=spec['rows'] or list(n.columns)!=spec['sample_columns']:raise ValueError('Candidate dimensions differ')
        if not np.isfinite(x).all() or (x<0).any() or not np.isclose(x.max(),spec['maximum']):raise ValueError('Candidate value audit differs')
        if not np.isclose(n.sum().min(),spec['minimum_column_sum']) or not np.isclose(n.sum().max(),spec['maximum_column_sum']):raise ValueError('Candidate column sums differ')
    if candidate['statistical_analysis_executed'] or candidate['primary_votes_added']!=0 or candidate['juvenile_standard_TPM_consistency_passed']:
        raise ValueError('Unresolved candidate promoted to analysis')
    return dict(queries=len(summary['queries']),unique_PMIDs=len(union),retrieved_abstract_records=len(seen),
        successful_payloads=payloads,recorded_failures=failures,quantitative_candidate_gate=candidate['quantitative_use'])


if __name__=='__main__':print(json.dumps(audit(ROOT),indent=2))
