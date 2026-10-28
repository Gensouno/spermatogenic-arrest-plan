#!/usr/bin/env python3
"""Reconcile independent-context tables and optionally reproduce in a fresh folder."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import pandas as pd
from run_independent_contexts import ROOT, COHORTS, digest

TABLES=['samples.csv','input_qc.csv','full_statistics.csv.gz','priority_effects.csv.gz',
        'fit_summary.csv','PCA.csv','gene_annotation.csv.gz','normalization.csv',
        'runtime_R_versions.csv','priority_stability.csv','priority_summary.csv',
        'discovery_inventory.csv','eligibility_audit.csv']
JSON_FILES=['summary.json','source_audit.json','discovery_summary.json','source_metadata.json']
FIGURES=['priority_context.png','sample_PCA.png']


def check(out):
    if len(pd.read_csv(out/'discovery_inventory.csv'))!=102 or len(pd.read_csv(out/'eligibility_audit.csv'))!=9:raise ValueError('Discovery review denominators')
    if len(json.loads((out/'source_metadata.json').read_text()))!=9:raise ValueError('Selected metadata evidence')
    samples=pd.read_csv(out/'samples.csv');fits=pd.read_csv(out/'fit_summary.csv')
    if len(samples)!=24 or samples.duplicated(['cohort','matrix_sample']).any():raise ValueError('Sample identity/count')
    for cohort in COHORTS:
        s=samples[samples.cohort==cohort];f=fits[fits.cohort==cohort]
        required={'full'}|{'omit_'+v for v in s.matrix_sample}
        if len(s)!=6 or s.groupby('group').size().to_dict()!={'Control':3,'Perturbation':3}:raise ValueError('Control/perturbation design')
        if len(f)!=7 or set(f.variant)!=required:raise ValueError('Omission suite incomplete')
    stats=pd.read_csv(out/'full_statistics.csv.gz');effects=pd.read_csv(out/'priority_effects.csv.gz')
    if stats.duplicated(['cohort','gene']).any() or effects.duplicated(['cohort','variant','symbol']).any():raise ValueError('Duplicate test or candidate')
    for c,d in stats.groupby('cohort'):
        if len(d)!=int(fits[(fits.cohort==c)&(fits.variant=='full')].tested_genes.iloc[0]):raise ValueError('Test family denominator')
        for field in ['P.Value','BH','BY']:
            if not np.isfinite(d[field]).all() or not d[field].between(0,1).all():raise ValueError('Invalid probability')
        p=d['P.Value'].to_numpy();order=np.argsort(p);n=len(p)
        bh=np.empty(n);bh[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1])
        by=np.minimum(1,bh*np.sum(1/np.arange(1,n+1)))
        if not np.allclose(d.BH,bh,rtol=1e-9,atol=1e-12) or not np.allclose(d.BY,by,rtol=1e-9,atol=1e-12):raise ValueError('Complete-family adjustment differs')
    stable=pd.read_csv(out/'priority_stability.csv')
    if len(stable)!=144 or stable.duplicated(['cohort','symbol']).any():raise ValueError('Fixed list denominator')
    for row in stable.itertuples():
        d=effects[(effects.cohort==row.cohort)&(effects.symbol==row.symbol)]
        s=samples[samples.cohort==row.cohort]
        complete=set(d.variant)==({'full'}|{'omit_'+v for v in s.matrix_sample})
        if row.complete!=complete or row.variants_measured!=len(d):raise ValueError('Missingness reconciliation')
        if row.all_negative!=bool(complete and (d.logFC<0).all()):raise ValueError('Sign retention')
        if row.all_down_1p5fold!=bool(complete and (d.logFC<=-np.log2(1.5)).all()):raise ValueError('Magnitude retention')
    summary=json.loads((out/'summary.json').read_text())
    if summary['primary_added_votes']!=0:raise ValueError('Exploratory branch changed primary family')
    return dict(tables=len(TABLES),json_files=len(JSON_FILES),figures=len(FIGURES),fits=len(fits),fixed_rows=len(stable))


def reproduce(out):
    with tempfile.TemporaryDirectory(prefix='independent-contexts-') as folder:
        fresh=Path(folder)
        subprocess.run(['bash','envs/python.sh','scripts/run_independent_contexts.py','--output',str(fresh)],cwd=ROOT,check=True)
        subprocess.run(['bash','envs/python.sh','scripts/plot_independent_contexts.py','--output',str(fresh)],cwd=ROOT,check=True)
        subprocess.run(['bash','envs/python.sh','scripts/audit_independent_metadata.py','--output',str(fresh)],cwd=ROOT,check=True)
        checked=[]
        for name in TABLES:
            if name in ['discovery_inventory.csv','eligibility_audit.csv']:continue
            a=pd.read_csv(out/name);b=pd.read_csv(fresh/name)
            pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=1e-9,atol=1e-10)
            checked.append(name)
        for name in ['summary.json','source_metadata.json']:
            if json.loads((out/name).read_text())!=json.loads((fresh/name).read_text()):raise ValueError('JSON differs: '+name)
        for name in COHORTS:
            pd.testing.assert_frame_equal(pd.read_csv(out/(name+'.tsv.gz'),sep='\t'),pd.read_csv(fresh/(name+'.tsv.gz'),sep='\t'))
        for name in FIGURES:
            if digest(out/name)!=digest(fresh/name):raise ValueError('Figure differs: '+name)
    record=dict(status='passed',scope='Fresh-directory numerical reconstruction from the same cached inputs; not independent animals or biological validation.',tables=checked,json_files=['summary.json','source_metadata.json'],figures=FIGURES,matrices=COHORTS,rtol=1e-9,atol=1e-10)
    (out/'reproduction.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reproduce',action='store_true');a=p.parse_args()
    out=ROOT/'results/independent_contexts';print(check(out))
    if a.reproduce:print(reproduce(out))
