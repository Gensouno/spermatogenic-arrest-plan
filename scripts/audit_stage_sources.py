#!/usr/bin/env python3
"""Reconcile cached Rif1 source records without extrapolating missing metadata."""
import argparse
import gzip
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from run_stage_specificity import ROOT
from audit_independent_metadata import parse

FIELDS={'strain','age','sex','tissue','replicate'}


def ena_records(payload):
    rows={}
    for sample in ET.fromstring(payload).findall('.//SAMPLE'):
        accession=sample.attrib.get('accession')
        attrs={a.findtext('TAG'):a.findtext('VALUE') for a in sample.findall('.//SAMPLE_ATTRIBUTE') if a.findtext('TAG') in FIELDS}
        if accession in rows:raise ValueError('Duplicate ENA sample')
        rows[accession]=dict(accession=accession,title=sample.findtext('TITLE'),attributes=attrs)
    return rows


def sra_records(payload):
    rows={}
    for package in ET.fromstring(payload).findall('.//EXPERIMENT_PACKAGE'):
        sample=package.find('SAMPLE');experiment=package.find('EXPERIMENT')
        ids=[e.text for e in sample.findall('./IDENTIFIERS/EXTERNAL_ID') if e.attrib.get('namespace')=='BioSample']
        if len(ids)!=1:raise ValueError('Expected one BioSample identity per SRA sample')
        accession=ids[0]
        if accession in rows:raise ValueError('Duplicate SRA BioSample')
        attrs={e.findtext('TAG'):e.findtext('VALUE') for e in sample.findall('.//SAMPLE_ATTRIBUTE') if e.findtext('TAG') in FIELDS}
        rows[accession]=dict(accession=accession,title=sample.findtext('TITLE'),attributes=attrs,
            experiment=experiment.attrib['accession'],sample=sample.attrib['accession'],
            library_name=experiment.findtext('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME'),
            runs=[e.attrib['accession'] for e in package.findall('./RUN_SET/RUN')])
    return rows


def main(out):
    raw=ROOT/'data/raw/stage_specificity'
    text=gzip.decompress((raw/'GSE294151_family.soft.gz').read_bytes()).decode()
    series,samples=parse(text)
    expected={};expected_experiments={};gsm=None
    for line in text.splitlines():
        if line.startswith('^SAMPLE = '):gsm=line.split(' = ',1)[1]
        if line.startswith('!Sample_relation = BioSample: '):
            match=re.search(r'SAMN\d+$',line)
            if not match:raise ValueError('Unexpected BioSample relation')
            expected[gsm]=match.group()
        if line.startswith('!Sample_relation = SRA: '):
            match=re.search(r'SRX\d+$',line)
            if not match:raise ValueError('Unexpected experiment relation')
            expected_experiments[gsm]=match.group()
    if len(expected)!=6 or len(set(expected.values()))!=6:raise ValueError('Expected six explicit GEO relations')
    batch=ena_records((raw/'Rif1_samples.xml').read_bytes())
    if not set(batch)<=set(expected.values()):raise ValueError('Unexpected batch accession')
    retrieval=json.loads((ROOT/'data/stage_source_retrieval.json').read_text())['individual']
    by_accession={a['accession']:a for a in retrieval}
    sra=sra_records((raw/'Rif1_SRA.xml').read_bytes())
    count=int(json.loads((raw/'SRA_count.json').read_text())['esearchresult']['count'])
    ids=json.loads((raw/'SRA_ids.json').read_text())['esearchresult']['idlist']
    if count!=6 or len(ids)!=6 or len(set(ids))!=6 or set(sra)!=set(expected.values()):raise ValueError('SRA count/accession reconciliation failed')
    rows=[]
    for gsm,accession in expected.items():
        a=by_accession[accession];p=raw/(accession+'.xml')
        row=dict(gsm=gsm,BioSample=accession,GEO_title=samples[gsm]['Sample_title'][0],ENA_status='not_retrieved',age=None)
        if 'error' in a:row['retrieval_error']=a['error']
        else:
            record=ena_records(p.read_bytes())
            if set(record)!={accession}:raise ValueError('Individual response differs from requested accession')
            e=record[accession];row.update(ENA_status='retrieved',ENA_title=e['title'],ENA_attributes=e['attributes'],ENA_age=e['attributes'].get('age'))
            if e['title']!=row['GEO_title']:raise ValueError('GEO and ENA titles disagree')
        s=sra[accession]
        if s['experiment']!=expected_experiments[gsm] or s['title']!=row['GEO_title']:raise ValueError('SRA/GEO identity differs')
        if row.get('ENA_age') and row['ENA_age']!=s['attributes'].get('age'):raise ValueError('ENA/SRA age disagreement')
        row.update(SRA_record=s,age=s['attributes'].get('age'))
        rows.append(row)
    summary=dict(recorded_on='2026-10-30',accession='GSE294151',expected_samples=len(expected),batch_returned=len(batch),
        individually_retrieved=sum(r['ENA_status']=='retrieved' for r in rows),
        SRA_records_reconciled=len(sra),ENA_age_p14_confirmed=sum(str(r.get('ENA_age')).lower()=='p14' for r in rows),
        age_p14_confirmed=sum(str(r['age']).lower()=='p14' for r in rows),
        age_unresolved=sum(r['age'] is None for r in rows),
        linked_GEO_citation=series.get('Series_pubmed_id',[]),
        samples=rows,limitations=['ENA returned two of six samples; four individual ENA404 records were resolved using the six GEO-linked SRA experiment records.','All six SRA records annotate age P14; archive accessions and biological-replicate labels do not independently establish distinct animals.','Exact Cre/allele, animal/pooling/batch structure, matrix-to-library crosswalk and phenotype remain unresolved.','GEO has no linked publication. Targeted literature search did not resolve a source paper, not proof that none exists.'])
    (out/'source_refresh_summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',default=str(ROOT/'results/stage_specificity'));a=parser.parse_args();main(Path(a.output))
