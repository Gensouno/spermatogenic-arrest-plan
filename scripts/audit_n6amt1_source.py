#!/usr/bin/env python3
"""Audit exact primary-cohort identities and preserve unresolved sampling fields."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import xml.etree.ElementTree as ET
from audit_independent_metadata import parse

ROOT=Path(__file__).resolve().parents[1]
FIELDS={'age','genotype','strain','tissue','sex','replicate'}


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def source_attributes(sample):
    return {e.findtext('TAG'):e.findtext('VALUE') for e in sample.findall('.//SAMPLE_ATTRIBUTE') if e.findtext('TAG') in FIELDS}


def require_identity(gsm,expected,experiment,sample):
    bio=[e.text for e in sample.findall('./IDENTIFIERS/EXTERNAL_ID') if e.attrib.get('namespace')=='BioSample']
    if len(bio)!=1 or bio[0]!=expected[gsm]['BioSample'] or experiment.attrib['accession']!=expected[gsm]['experiment']:
        raise ValueError('GEO/SRA identity mismatch')


def audit(root,out):
    out.mkdir(parents=True,exist_ok=True)
    for a in json.loads((root/'data/n6amt1_source_inputs.json').read_text())['assets']:
        p=root/a['path']
        if p.stat().st_size!=a['bytes'] or digest(p)!=a['sha256']:raise ValueError('Source identity differs: '+a['path'])
    raw=root/'data/raw/n6amt1_source';text=gzip.decompress((raw/'GSE283904_family.soft.gz').read_bytes()).decode()
    series,samples=parse(text);expected={};gsm=None
    for line in text.splitlines():
        if line.startswith('^SAMPLE = '):gsm=line.split(' = ',1)[1];expected[gsm]={}
        if line.startswith('!Sample_relation = BioSample: '):expected[gsm]['BioSample']=re.search(r'SAMN\d+$',line).group()
        if line.startswith('!Sample_relation = SRA: '):expected[gsm]['experiment']=re.search(r'SRX\d+$',line).group()
    if len(expected)!=8 or len({a['BioSample'] for a in expected.values()})!=8:raise ValueError('Expected eight GEO samples')
    source=ET.parse(raw/'N6amt1_SRA.xml');packages=source.findall('.//EXPERIMENT_PACKAGE')
    count=int(json.loads((raw/'SRA_count.json').read_text())['esearchresult']['count']);uids=json.loads((raw/'SRA_ids.json').read_text())['esearchresult']['idlist']
    if count!=8 or len(set(uids))!=8 or len(packages)!=8:raise ValueError('Incomplete SRA retrieval')
    old={a['sample_id']:a for a in csv.DictReader((root/'data/samples.tsv').open(),delimiter='\t')}
    rows=[];seen=set()
    for package in packages:
        e=package.find('EXPERIMENT');s=package.find('SAMPLE');matches=[k for k,v in expected.items() if v['experiment']==e.attrib['accession']]
        if len(matches)!=1:raise ValueError('Ambiguous experiment match')
        gsm=matches[0];require_identity(gsm,expected,e,s)
        if gsm in seen:raise ValueError('Repeated sample')
        seen.add(gsm);attrs=source_attributes(s);title=samples[gsm]['Sample_title'][0]
        if title!=s.findtext('TITLE'):raise ValueError('Source title disagreement')
        characteristics=samples[gsm].get('Sample_characteristics_ch1',[])
        age=[v.split(':',1)[1].strip() for v in characteristics if v.lower().startswith('age:')]
        rows.append(dict(gsm=gsm,title=title,BioSample=expected[gsm]['BioSample'],experiment=e.attrib['accession'],
            source_sample=s.attrib['accession'],GEO_genotype=';'.join(v.split(':',1)[1].strip() for v in characteristics if v.lower().startswith('genotype:')),
            SRA_genotype=attrs.get('genotype',''),GEO_age=';'.join(age),SRA_age=attrs.get('age',''),
            old_curation_age=old[gsm].get('age',''),source_age_status='not reported' if not age and not attrs.get('age') else 'reported',
            analysis_label='WT' if title.startswith('CTL-') else 'KO',literal_WT_genotype_verified=attrs.get('genotype','').lower() in ['wt','wild type','wild-type']))
    rows=sorted(rows,key=lambda a:a['gsm'])
    with (out/'sample_audit.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    pub=ET.parse(raw/'pubmed.xml');pmid=pub.findtext('.//PMID');crossref=json.loads((root/'paper/literature/Luo2025.json').read_text())['message']
    if pmid!='40079145' or crossref['DOI'].lower()!='10.1093/biolre/ioaf052':raise ValueError('Related article identity mismatch')
    report=dict(recorded_on='2026-10-31',accession='GSE283904',samples=8,controls=sum(a['title'].startswith('CTL-') for a in rows),perturbations=sum(a['title'].startswith('N6A-SKO-') for a in rows),
        source_age_missing=sum(a['source_age_status']=='not reported' for a in rows),literal_WT_controls_verified=sum(a['literal_WT_genotype_verified'] for a in rows),
        related_article=dict(PMID=pmid,DOI=crossref['DOI'],publication_status='published',role='Related biological study; exact GSE283904 linkage not verified from accessible article content',fulltext_retrieved=False),
        GEO_linked_citation=series.get('Series_pubmed_id',[]),primary_statistics_changed=False,
        limitations=['Old adult curation entries are not confirmed source ages.','CTL is a source control label; WT is an analysis alias, not a verified literal control genotype.','Exact Cre/alleles, animal/pooling/batch and sample-matched phenotype remain unresolved.','The related article abstract reports meiotic abnormalities but does not establish the accession-to-article link; publisher retrieval returned403.'])
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def freeze(root,out):
    target=root/'n6amt1_source_snapshot'
    if (target/'manifest.json').exists():raise ValueError('Snapshot frozen; choose a new name')
    target.mkdir(exist_ok=True)
    for name in ['sample_audit.csv','summary.json']:shutil.copy2(out/name,target/name)
    sources=['scripts/audit_n6amt1_source.py','scripts/audit_independent_metadata.py','scripts/run_independent_contexts.py','tests/test_n6amt1_source.py']
    specs=['data/n6amt1_source_contract.json','data/n6amt1_source_inputs.json','docs/N6AMT1_SOURCE.md']
    manifest=dict(scope='Current source identity and unresolved fields; not biological validation or new statistics.',source_sha256={p:digest(root/p) for p in sources},input_specification_sha256={p:digest(root/p) for p in specs},files=[dict(path=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(target.iterdir()) if p.is_file()])
    (target/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--freeze',action='store_true');p.add_argument('--output',default=str(ROOT/'results/n6amt1_source'));a=p.parse_args();out=Path(a.output);print(audit(ROOT,out))
    if a.freeze:freeze(ROOT,out)
