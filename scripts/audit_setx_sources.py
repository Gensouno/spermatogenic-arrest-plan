#!/usr/bin/env python3
"""Reconcile frozen Setx source identities; never infer animals from accessions."""
import gzip
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def check_inputs(root):
    for item in json.loads((root/'data/setx_source_audit_inputs.json').read_text())['assets']:
        p=Path(item['path']);target=root/p
        if p.is_absolute() or '..' in p.parts or not target.resolve().is_relative_to(root.resolve()):
            raise ValueError('Unsafe source input')
        if target.stat().st_size!=item['bytes'] or digest(target)!=item['sha256']:
            raise ValueError('Source identity differs: '+str(p))


def parse_soft(text):
    series={};samples={};current=None
    for line in text.splitlines():
        if line.startswith('^'):
            current=None
            if line.startswith('^SAMPLE = '):
                gsm=line.split(' = ',1)[1]
                if gsm in samples:raise ValueError('Duplicate GEO sample record')
                samples[gsm]={};current=samples[gsm]
        if ' = ' not in line:continue
        k,v=line.split(' = ',1)
        if k.startswith('!Series_'):series.setdefault(k,[]).append(v)
        if current is not None and k.startswith('!Sample_'):current.setdefault(k,[]).append(v)
    return series,samples


def exact_relation(sample,prefix,pattern):
    values=[v for v in sample.get('!Sample_relation',[]) if v.startswith(prefix+':')]
    ids=[re.findall(pattern,v) for v in values]
    if len(ids)!=1 or len(ids[0])!=1:raise ValueError('Ambiguous '+prefix+' relation')
    return ids[0][0]


def reconcile(frozen,current,ena,expected_count,design):
    gsms=set(frozen)
    if set(current)!=gsms or set(design.gsm)!=gsms or design.gsm.duplicated().any():
        raise ValueError('GEO/design sample set differs')
    if len(ena)!=expected_count or set(ena.sample_alias)!=gsms:
        raise ValueError('ENA count/sample reconciliation failed')
    for field in ['run_accession','experiment_accession','sample_accession','sample_alias']:
        if ena[field].isna().any() or ena[field].duplicated().any():raise ValueError('Missing or duplicate ENA '+field)
    ena=ena.set_index('sample_alias');design=design.set_index('gsm');rows=[]
    for gsm in sorted(gsms):
        old=frozen[gsm];new=current[gsm];e=ena.loc[gsm];d=design.loc[gsm]
        for key in ['!Sample_title','!Sample_characteristics_ch1','!Sample_supplementary_file_1','!Sample_supplementary_file_2']:
            if old.get(key)!=new.get(key):raise ValueError('GEO source label/file changed: '+gsm+' '+key)
        title=new['!Sample_title'][0]
        if e.sample_title!=title:raise ValueError('ENA/GEO title mismatch')
        match=re.fullmatch(r'(wt|SetxKO) testis, (10|12|14|16|18)dpp, rep([1-4])',title)
        if not match:raise ValueError('Unexpected source title')
        age='P'+match[2];group='Control' if match[1]=='wt' else 'KO'
        if age!=d.age or group!=d.group:raise ValueError('Source/design label mismatch')
        biosample=exact_relation(new,'BioSample',r'SAMN\d+')
        experiment=exact_relation(new,'SRA',r'SRX\d+')
        if biosample!=e.sample_accession or experiment!=e.experiment_accession:
            raise ValueError('GEO/ENA accession mismatch')
        if str(e.tax_id)!='10090' or e.library_strategy!='RNA-Seq':raise ValueError('Wrong organism or assay')
        rows.append(dict(gsm=gsm,age=age,source_group_label=group,source_genotype=match[1],replicate_label=int(match[3]),
            filename_source_code=str(d.source_code),barcode_token=d.barcode_token,
            biosample=biosample,sra_experiment=experiment,sra_sample=e.secondary_sample_accession,
            run=e.run_accession,library_name=e.library_name,layout=e.library_layout,instrument=e.instrument_model,
            source_title=title,read_count=int(e.read_count),base_count=int(e.base_count),
            source_identity_reconciled=True,individual_animal_identity_confirmed=False))
    return pd.DataFrame(rows)


class SourceText(HTMLParser):
    def __init__(self):super().__init__();self.parts=[];self.skip=0;self.metadata={}
    def handle_starttag(self,tag,attrs):
        if tag in ['script','style']:self.skip+=1
        a=dict(attrs)
        if tag=='meta' and a.get('name'):self.metadata.setdefault(a['name'],[]).append(a.get('content',''))
    def handle_endtag(self,tag):
        if tag in ['script','style']:self.skip-=1
    def handle_data(self,data):
        if not self.skip:self.parts.append(data)


def publication_evidence(raw):
    p=SourceText();p.feed((raw/'preprint.html').read_text());text=re.sub(r'\s+',' ',' '.join(p.parts))
    # Verify bounded source statements before assigning their scoped annotations.
    required=['GSE224875','18 dpp (four biological replicates)','induced in a C57BL/6J',
              'C3H/HeJ','30 generations','100-bp pair-end reads']
    if any(v not in text for v in required):raise ValueError('Expected preprint evidence changed')
    pub=ET.parse(raw/'pubmed.xml').getroot()
    if pub.findtext('.//PMID')!='37801077':raise ValueError('PubMed identity differs')
    abstract=' '.join(''.join(e.itertext()) for e in pub.findall('.//AbstractText'))
    if 'absence of spermatids and mature spermatozoa' not in abstract or 'mid-prophase' not in abstract:
        raise ValueError('Published abstract phenotype evidence changed')
    doi=[e.text for e in pub.findall('.//ArticleId') if e.get('IdType')=='doi']
    if doi!=['10.1530/REP-23-0166']:raise ValueError('Published DOI mismatch')
    api=json.loads((raw/'preprint_metadata_doi.json').read_text())['collection']
    if len(api)!=1 or api[0]['doi']!='10.1101/2023.04.12.536672' or api[0]['version']!='1':
        raise ValueError('Preprint identity/version changed')
    if p.metadata.get('citation_id')!=['2023.04.12.536672v1']:raise ValueError('HTML source is not the specified preprint version')
    jats=ET.parse(raw/'preprint_source.xml').getroot()
    names=[e.findtext('given-names')+' '+e.findtext('surname') for e in jats.findall('.//contrib[@contrib-type="author"]/name')]
    if p.metadata.get('citation_author')!=names or len(names)!=8:raise ValueError('Versioned HTML/JATS author identity differs')
    return dict(published_PMID='37801077',published_DOI=doi[0],published_year=2023,
        preprint_DOI=api[0]['doi'],preprint_version=1,preprint_document_authors=names,
        mutant_allele='Setx spcar3; ENU-induced mutant allele, not evidence of an engineered knockout',
        allele_association='GEO series summary names spcar3; the openly accessible preprint explicitly names GSE224875.',
        histology_scope='Published abstract reports absence of spermatids/mature sperm and mid-prophase meiotic arrest. Author-reported allele-level evidence; not matched histology for each reanalyzed library.',
        preprint_RNAseq_age_dpp=18,preprint_reported_biological_replicates=4,
        other_age_independence='The P18 statement is not extended to P10-P16. No animal manifest or pooling map was established.',
        preprint_background='C3H/HeJ congenic line after more than 30 backcross generations; reported study context, not an individual-library strain certificate.',
        published_fulltext_locally_retrieved=False,
        access_limit='Publisher HTML/PDF could not be downloaded in this audit; published phenotype was checked in PubMed. Detailed methods checked in preprint v1, explicitly not equated with the final text.',
        preprint_latest_metadata_caveat='The bioRxiv API response says published=NA and omits Sabrina Petri; the version-1 HTML/JATS and Crossref include eight authors. Use the versioned document and Crossref for author identity, not the inconsistent API summary. Independent PubMed/Crossref identify a published article; publication versions remain distinct.')


def main():
    check_inputs(ROOT);raw=ROOT/'data/raw/setx_source_audit';out=ROOT/'results/setx_source_audit';out.mkdir(parents=True,exist_ok=True)
    _,old=parse_soft(gzip.decompress((ROOT/'data/raw/replicability_discovery/GSE224875_family.soft.gz').read_bytes()).decode())
    series,current=parse_soft((raw/'geo_series.soft').read_text())
    ena=pd.read_csv(raw/'ena_runs.tsv',sep='\t',dtype=str).fillna('')
    count=int((raw/'ena_count.txt').read_text().splitlines()[-1])
    design=pd.read_csv(ROOT/'setx_interaction_snapshot/sample_design.csv',dtype=str)
    table=reconcile(old,current,ena,count,design)
    table.to_csv(out/'sample_identity.csv',index=False)
    pair=table[table.filename_source_code.duplicated(False)]
    pair.to_csv(out/'ambiguous_pair.csv',index=False)
    # Check exact nested links for both libraries; record attributes, not inferred animals.
    details=[]
    for gsm in ['GSM7034763','GSM7034783']:
        row=table[table.gsm==gsm].iloc[0]
        sr=ET.parse(raw/(gsm+'_sample.xml')).getroot().find('SAMPLE')
        rr=ET.parse(raw/(gsm+'_run.xml')).getroot().find('RUN')
        er=ET.parse(raw/('experiment_WT_P16.xml' if row.source_group_label=='Control' else 'experiment_Mut_P16.xml')).getroot().find('EXPERIMENT')
        if sr.get('accession')!=row.biosample or rr.get('accession')!=row.run or er.get('accession')!=row.sra_experiment:
            raise ValueError('Ambiguous-pair XML identity mismatch')
        if rr.find('EXPERIMENT_REF').get('accession')!=row.sra_experiment:raise ValueError('Run/experiment link mismatch')
        attrs={e.findtext('TAG'):e.findtext('VALUE') for e in sr.findall('.//SAMPLE_ATTRIBUTE')}
        details.append(dict(gsm=gsm,sample_attributes=attrs,animal_identity_confirmed=False,
            interpretation='Distinct archive records; repeated filename code remains unresolved. No label correction inferred.'))
    (out/'pair_attribute_audit.json').write_text(json.dumps(details,indent=2)+'\n')
    publication=publication_evidence(raw)
    if not any('spcar3' in v for v in series.get('!Series_summary',[])):raise ValueError('GEO allele link missing')
    (out/'publication_evidence.json').write_text(json.dumps(publication,indent=2)+'\n')
    summary=dict(scope='Primary source identity, author-reported phenotype and version-aware publication audit; not biological validation.',
        frozen_GEO_samples=len(old),current_GEO_samples=len(current),ENA_count=count,ENA_runs=len(ena),reconciled_samples=len(table),
        unique_biosamples=int(table.biosample.nunique()),unique_experiments=int(table.sra_experiment.nunique()),unique_runs=int(table.run.nunique()),
        paired_libraries=int((table.layout=='PAIRED').sum()),HiSeq2500_libraries=int((table.instrument=='Illumina HiSeq 2500').sum()),
        repeated_filename_code_samples=len(pair),independent_animals_confirmed=False,
        source_labels_unchanged=True,barcode_confounding_resolved=False,
        conclusion='All sample/run identities reconcile; separate accessions do not resolve filename code 7844 or prove independent animals. The original RNA-seq source reports biological replication at P18 in preprint v1, without a complete cross-age animal/pooling manifest. Source KO label is retained as a mutant-contrast alias; the biological allele is spcar3.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
