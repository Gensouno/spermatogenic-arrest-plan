#!/usr/bin/env python3
"""Check the complete competitive family, coverage, source fits and reproduction."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import pandas as pd
from run_competitive_set import ROOT, COHORTS, digest

TABLES=['camera_results.csv','gene_universe.csv.gz','priority_coverage.csv',
        'fit_reconciliation.csv','runtime_R_versions.csv']
JSON_FILES=['summary.json','runtime_Python_versions.json']
TEXT_FILES=['installed_camera_source.txt']
FIGURES=['correlation_sensitivity.png']


def adjust_family(frame):
    """Independently reconcile the eight-test family, never separately by policy."""
    expected={(c,p) for c in COHORTS for p in ['estimated','fixed_0p01']}
    if len(frame)!=8 or set(zip(frame.cohort,frame.policy))!=expected:
        raise ValueError('Complete eight-test family required')
    p=frame.PValue.to_numpy(float)
    if not np.isfinite(p).all() or ((p<0)|(p>1)).any():raise ValueError('Invalid P values')
    order=np.argsort(p,kind='stable');n=len(p);bh=np.empty(n)
    bh[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1])
    by=np.minimum(1,bh*np.sum(1/np.arange(1,n+1)))
    return bh,by


def validate_results(results):
    bh,by=adjust_family(results)
    if not np.allclose(results.BH,bh,rtol=1e-9,atol=1e-30) or not np.allclose(results.BY,by,rtol=1e-9,atol=1e-30):
        raise ValueError('Joint-family adjustment differs')
    if not results.fixed_priorities.eq(36).all():raise ValueError('Fixed priority denominator changed')
    if not results.Direction.isin(['Up','Down']).all():raise ValueError('Direction invalid')
    reverse=results.Direction.map({'Up':'Down','Down':'Up'})
    if not reverse.eq(results.reverse_Direction).all() or not np.allclose(results.PValue,results.reverse_PValue,rtol=1e-9,atol=1e-30):
        raise ValueError('Contrast reversal failed')
    for row in results.itertuples():
        if not np.isfinite(row.correlation_raw_or_assumed) or not -1<=row.correlation_raw_or_assumed<=1:
            raise ValueError('Correlation invalid')
        expected_df=4 if row.policy=='estimated' else row.universe_genes-2
        if row.camera_df!=expected_df:raise ValueError('Correlation policy degrees of freedom differ')
        if row.policy=='fixed_0p01' and not np.isclose(row.correlation_raw_or_assumed,.01):raise ValueError('Fixed correlation differs')
        vif=max(1,1+(row.NGenes-1)*row.correlation_raw_or_assumed)
        if not np.isfinite(vif) or not np.isclose(row.variance_inflation_used,vif):raise ValueError('Variance inflation differs')
        if row.background_genes!=row.universe_genes-row.NGenes:raise ValueError('Competitive complement differs')
    for _,d in results.groupby('cohort'):
        for field in ['NGenes','universe_genes','background_genes','priority_mean_logFC','background_mean_logFC']:
            if d[field].nunique()!=1:raise ValueError('Policy-dependent universe or gene effects')


def check(out):
    results=pd.read_csv(out/'camera_results.csv');validate_results(results)
    universe=pd.read_csv(out/'gene_universe.csv.gz');coverage=pd.read_csv(out/'priority_coverage.csv')
    audit=pd.read_csv(out/'fit_reconciliation.csv')
    priorities=set(pd.read_csv(ROOT/'replicability_snapshot/stable_descriptive_priorities.csv').gene)
    primary=set(pd.read_csv(ROOT/'replicability_snapshot/all_primary_gene_stage_context.csv.gz').gene)
    original=pd.read_csv(ROOT/'independent_contexts_snapshot/full_statistics.csv.gz')
    if universe.duplicated(['cohort','gene']).any() or universe.duplicated(['cohort','symbol']).any():raise ValueError('Duplicate universe identity')
    if len(coverage)!=144 or coverage.duplicated(['cohort','symbol']).any():raise ValueError('Coverage inventory differs')
    if len(audit)!=4 or set(audit.cohort)!=set(COHORTS):raise ValueError('Fit inventory differs')
    for cohort in COHORTS:
        u=universe[universe.cohort==cohort];cov=coverage[coverage.cohort==cohort]
        previous=original[(original.cohort==cohort)&original.globally_unique_symbol.eq(True)&original.symbol.isin(primary)]
        if set(u.gene)!=set(previous.gene) or set(u.symbol)!=set(previous.symbol):raise ValueError('Declared universe differs')
        if not u.priority.eq(u.symbol.isin(priorities)).all() or set(cov.symbol)!=priorities:raise ValueError('Priority membership differs')
        if not cov.available.eq(cov.symbol.isin(u.symbol)).all():raise ValueError('Missingness differs')
        row=results[results.cohort==cohort].iloc[0];a=audit[audit.cohort==cohort].iloc[0]
        if row.NGenes!=int(cov.available.sum()) or row.universe_genes!=len(u):raise ValueError('Coverage denominator differs')
        if a.source_tested_genes!=len(original[original.cohort==cohort]) or a.measured_priorities!=row.NGenes or a.universe_genes!=len(u):raise ValueError('Fit reconciliation denominator differs')
        if (a[['logFC_max_error','AveExpr_max_error','t_max_error','PValue_max_error']].astype(float)>1e-8).any():raise ValueError('Old full fits not reproduced')
        if (a.libraries,a.controls,a.perturbations,a.residual_df)!=(6,3,3,4):raise ValueError('Library design differs')
        if not np.isclose(row.priority_mean_logFC,u.loc[u.priority,'logFC'].mean()) or not np.isclose(row.background_mean_logFC,u.loc[~u.priority,'logFC'].mean()):raise ValueError('Effect summary differs')
    summary=json.loads((out/'summary.json').read_text())
    if summary['family_tests']!=8 or summary['primary_discoveries_changed'] is not False:raise ValueError('Primary family boundary differs')
    for policy,d in results.groupby('policy'):
        expected=dict(contexts=len(d),direction_down=int(d.Direction.eq('Down').sum()),
            BY_below_0p05=d.loc[d.BY.lt(.05),'cohort'].tolist(),BH_below_0p05=d.loc[d.BH.lt(.05),'cohort'].tolist())
        if summary['results'][policy]!=expected:raise ValueError('Summary differs')
    return dict(tests=len(results),coverage_records=len(coverage),universe_rows=len(universe),source_fits=len(audit))


def reproduce(out):
    with tempfile.TemporaryDirectory(prefix='competitive-set-') as folder:
        fresh=Path(folder)
        subprocess.run(['bash','envs/python.sh','scripts/run_competitive_set.py','--output',str(fresh)],cwd=ROOT,check=True)
        subprocess.run(['bash','envs/python.sh','scripts/plot_competitive_set.py','--output',str(fresh)],cwd=ROOT,check=True)
        check(fresh)
        for name in TABLES:
            pd.testing.assert_frame_equal(pd.read_csv(out/name),pd.read_csv(fresh/name),check_exact=False,rtol=1e-9,atol=1e-12)
        for name in JSON_FILES:
            if json.loads((out/name).read_text())!=json.loads((fresh/name).read_text()):raise ValueError('JSON reproduction differs')
        for name in TEXT_FILES+FIGURES:
            if digest(out/name)!=digest(fresh/name):raise ValueError('Source/figure reproduction differs: '+name)
        for cohort in COHORTS:
            name=cohort+'.tsv.gz'
            pd.testing.assert_frame_equal(pd.read_csv(out/'prepared'/name,sep='\t'),pd.read_csv(fresh/'prepared'/name,sep='\t'))
    record=dict(status='passed',scope='Fresh cached-source reconstruction and numerical rerun; not independent biological validation.',
        tables=TABLES,json_files=JSON_FILES,text_files=TEXT_FILES,figures=FIGURES,matrices=4,rtol=1e-9,atol=1e-12)
    (out/'reproduction.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--reproduce',action='store_true');a=parser.parse_args()
    out=ROOT/'results/competitive_set';print(check(out))
    if a.reproduce:print(reproduce(out))
