#!/usr/bin/env python3
"""Check annotation-vs-count columns and select genome RSEM tables, never piRNA."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
EEF_COLUMNS=[f'Ee_{group}{i}' for group in ['con','cKO'] for i in [1,2,3]]


def read_eef(path):
    d=pd.read_csv(path,sep='\t')
    if not {'gene_id','gene_name',*EEF_COLUMNS}<=set(d) or d.gene_id.isna().any() or d.gene_id.duplicated().any():
        raise ValueError('Eef1g count schema/identifiers differ')
    if [c for c in d if c.startswith('Ee_')]!=EEF_COLUMNS:raise ValueError('Unexpected count sample column/order')
    x=d[EEF_COLUMNS].to_numpy(dtype=float)
    if not np.isfinite(x).all() or (x<0).any() or not np.allclose(x,np.round(x),atol=1e-8,rtol=0):
        raise ValueError('Eef1g source is not nonnegative integer counts')
    return d


def select_genome_members(members,gsms):
    selected={};pirna=set()
    for m in members:
        if not m.isfile():raise ValueError('Unexpected non-file archive member')
        if '/' in m.name or '\\' in m.name:raise ValueError('Unexpected nested archive path')
        if '_mm10_piRNA.' in m.name:
            match=re.fullmatch(r'(GSM[0-9]+)_.+_mm10_piRNA\.trimmed\.genes\.results\.gz',m.name)
            if not match or match[1] not in gsms or match[1] in pirna:raise ValueError('Unexpected or duplicate piRNA member')
            pirna.add(match[1]);continue
        match=re.fullmatch(r'(GSM[0-9]+)_.+_GRCm38\.trimmed\.genes\.results\.gz',m.name)
        if not match or match[1] not in gsms or match[1] in selected:raise ValueError('Unexpected or duplicate genome member')
        selected[match[1]]=m
    if set(selected)!=set(gsms):raise ValueError('Genome/member GSM reconciliation differs')
    if pirna!=set(gsms):raise ValueError('piRNA exclusion inventory differs')
    return selected


def main():
    out=ROOT/'results/new_models';out.mkdir(parents=True,exist_ok=True)
    prepared=ROOT/'data/interim/new_models';prepared.mkdir(parents=True,exist_ok=True)
    meta=ROOT/'results/replicability_discovery'
    eef=read_eef(ROOT/'data/raw/replicability_discovery/GSE307441_raw_counts.txt.gz')
    eef.set_index('gene_id')[EEF_COLUMNS].to_csv(prepared/'Eef1g_counts.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
    eef[['gene_id','gene_name']].rename(columns={'gene_name':'symbol'}).to_csv(out/'Eef1g_source_annotation.csv',index=False)
    evidence=[]
    source=json.loads((meta/'GSE307441_metadata.json').read_text())
    for i,column in enumerate(EEF_COLUMNS):
        gsm=f'GSM{9224238+i}';s=source['samples'][gsm];group='Control' if i<3 else 'KO';rep=i%3+1
        expected=f'{"Control" if group=="Control" else "Eef1g cKO"} leptotene/zygotene spermatocytes_rep{rep}'
        if s['title']!=[expected] or s['organism_ch1']!=['Mus musculus']:raise ValueError('Eef1g source sample identity differs')
        evidence.append({'cohort':'Eef1g_LZ_P12','gse':'GSE307441','gsm':gsm,'matrix_sample':column,'group':group,
            'age':'P12','source_title':expected,'source_genotype':'Wild type' if group=='Control' else 'Eef1g cKO',
            'scope':'FACS L/Z; individual animal independence and stage proportions unconfirmed'})
    setx=json.loads((meta/'GSE224875_metadata.json').read_text())
    if len(setx['samples'])!=40:raise ValueError('Setx source sample count differs')
    inventory=[];matrices={k:{} for k in ['TPM','expected_count','effective_length']};genes=None
    with tarfile.open(ROOT/'data/raw/replicability_discovery/GSE224875_RAW.tar') as archive:
        members=archive.getmembers();selected=select_genome_members(members,setx['samples'])
        if len(members)!=80:raise ValueError('Expected 40 genome + 40 piRNA members')
        for gsm,m in sorted(selected.items()):
            compressed=archive.extractfile(m).read();d=pd.read_csv(io.BytesIO(gzip.decompress(compressed)),sep='\t')
            needed=['expected_count','effective_length','TPM']
            if not {'gene_id',*needed}<=set(d) or d.gene_id.isna().any() or d.gene_id.duplicated().any():raise ValueError('Invalid RSEM gene schema')
            values=d[needed].to_numpy(dtype=float)
            if not np.isfinite(values).all() or (values<0).any():raise ValueError('Invalid RSEM abundance')
            if genes is None:genes=d.gene_id
            if not genes.equals(d.gene_id):raise ValueError('RSEM gene ordering/identity differs')
            total=float(d.TPM.sum())
            if abs(total-1e6)>500:raise ValueError('Unexpected RSEM TPM total')
            for k in matrices:matrices[k][gsm]=d[k].to_numpy()
            s=setx['samples'][gsm];title=s['title'][0]
            match=re.fullmatch(r'(wt|SetxKO) testis, (10|12|14|16|18)dpp, rep([1-4])',title)
            if not match or s['organism_ch1']!=['Mus musculus']:raise ValueError('Setx source sample label differs')
            group='Control' if match[1]=='wt' else 'KO';age='P'+match[2]
            expected_genotype='genotype: WT' if group=='Control' else 'genotype: SetxKO'
            if expected_genotype not in s['characteristics_ch1']:raise ValueError('Setx genotype contradicts title')
            evidence.append({'cohort':'Setx_'+age,'gse':'GSE224875','gsm':gsm,'matrix_sample':gsm,'group':group,
                'age':age,'source_title':title,'source_genotype':expected_genotype.split(': ')[1],
                'scope':'Whole testis; RSEM estimates/TPM, animal independence unconfirmed'})
            inventory.append({'gsm':gsm,'member':m.name,'compressed_bytes':len(compressed),
                'member_sha256':hashlib.sha256(compressed).hexdigest(),'genes':len(d),'TPM_sum':total,
                'expected_count_sum':float(d.expected_count.sum()),
                'fraction_noninteger_expected_count':float((np.abs(d.expected_count-np.round(d.expected_count))>1e-6).mean())})
    for k,columns in matrices.items():
        frame=pd.DataFrame(columns,index=pd.Index(genes,name='gene'))
        if k=='TPM' and frame.T.duplicated().any():raise ValueError('Identical genome TPM sample vectors')
        frame.to_csv(prepared/f'Setx_{k}.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
    evidence=pd.DataFrame(evidence)
    sizes=evidence.groupby(['cohort','group']).size()
    if not ((sizes.loc['Eef1g_LZ_P12']==3).all() and (sizes.drop('Eef1g_LZ_P12',level=0)==4).all()):raise ValueError('New model group balance differs')
    evidence.to_csv(out/'sample_evidence.csv',index=False);pd.DataFrame(inventory).to_csv(out/'Setx_archive_audit.csv',index=False)
    report={'Eef1g_gene_rows':len(eef),'Eef1g_count_columns':EEF_COLUMNS,'Eef1g_samples':6,
        'Setx_genome_samples':40,'Setx_gene_rows':len(genes),'Setx_piRNA_members_excluded':40,
        'Setx_units':'Separate RSEM expected_count, effective_length and TPM. Descriptive analysis uses TPM only; estimates not rounded into raw counts.',
        'sample_identity':'GEO titles/genotypes reconciled; Setx genome filename GSM exact. Eef1g source control/cKO and replicate suffix to source title, no gene-based relabelling.',
        'limitations':['No source-independent confirmation of animals, pooling, or hidden batches.',
            'Setx sampling ages from source titles; original library/annotation details not reconstructed.',
            'Sorted L/Z mixture and whole testis are distinct contexts, never additional votes in the existing four-unit test.']}
    (out/'preparation_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
