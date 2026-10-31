#!/usr/bin/env python3
"""Audit eligible targets, calipers, equal-target weighting and cached-input reproduction."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import pandas as pd
from run_stage_specificity import ROOT,digest

TABLES=['gene_covariates.csv.gz','covariate_reconstruction.csv','target_eligibility.csv',
        'neighborhoods.csv.gz','paired_comparisons.csv','specificity_summary.csv','stage_descriptor_audit.csv',
        'resolution_eligibility.csv','resolution_neighborhoods.csv.gz','resolution_pairs.csv','resolution_summary.csv']
JSON_FILES=['design_summary.json','runtime_versions.json','source_refresh_summary.json','interpretation_summary.json','resolution_design.json']
FIGURES=['matched_backgrounds.png','stage_resolution.png']


def check(out):
    eligibility=pd.read_csv(out/'target_eligibility.csv');edges=pd.read_csv(out/'neighborhoods.csv.gz')
    pairs=pd.read_csv(out/'paired_comparisons.csv');summary=pd.read_csv(out/'specificity_summary.csv')
    keys=['cohort','reference','abundance_caliper','comparison']
    if len(eligibility)!=576 or eligibility.duplicated(keys[:3]+['symbol']).any():raise ValueError('Priority denominator/inventory differs')
    if len(summary)!=32 or summary.duplicated(keys).any():raise ValueError('Scenario inventory differs')
    if edges.duplicated(keys+['target','background']).any():raise ValueError('Duplicated neighborhood member')
    priorities=set(eligibility.symbol)
    if len(priorities)!=36 or edges.background.isin(priorities).any():raise ValueError('Priority leaked into background')
    if (edges.abundance_gap>edges.abundance_caliper+1e-11).any():raise ValueError('Abundance caliper violated')
    strict=edges[edges.comparison=='stage_abundance_dominance']
    if (strict.target_peak_group!=strict.background_peak_group).any() or (strict.dominance_gap>.2+1e-11).any():raise ValueError('Stage/dominance caliper violated')
    for key,d in pairs.groupby(keys,sort=True):
        s=summary
        for k,value in zip(keys,key):s=s[s[k]==value]
        if len(s)!=1:raise ValueError('Summary key ambiguity')
        row=s.iloc[0]
        if len(d)!=row.matched_priorities or not d.background_genes.between(5,20).all():raise ValueError('Matched denominator/neighbor count')
        values=[d.target_logFC.mean(),d.background_mean_logFC.mean(),(d.target_logFC-d.background_mean_logFC).mean(),
                (d.target_logFC-d.background_mean_logFC).median()]
        observed=[row.priority_mean_logFC,row.background_mean_logFC,row.mean_logFC_difference,row.median_pairwise_logFC_difference]
        if not np.allclose(values,observed,rtol=1e-10,atol=1e-11):raise ValueError('Equal-target effect reconciliation')
    for key,d in edges.groupby(keys+['target'],sort=True):
        p=pairs
        for k,value in zip(keys+['target'],key):p=p[p[k]==value]
        if len(p)!=1 or not np.isclose(d.background_logFC.mean(),p.background_mean_logFC.iloc[0],rtol=1e-10,atol=1e-11):raise ValueError('Neighborhood effect mismatch')
    for key,d in pairs.groupby(keys[:3],sort=True):
        targets=[set(x.target) for _,x in d.groupby('comparison')]
        if len(targets)!=2 or targets[0]!=targets[1]:raise ValueError('Different target sets between schemes')
    audit=pd.read_csv(out/'covariate_reconstruction.csv')
    if audit.voom_AveExpr_max_absolute_error.max()>1e-10 or audit.control_normalization_shift_spread.max()>1e-10:raise ValueError('Voom/control abundance audit')
    source=json.loads((out/'source_refresh_summary.json').read_text())
    if source['expected_samples']!=6 or source['SRA_records_reconciled']!=6 or source['age_p14_confirmed']!=6:raise ValueError('Rif1 source reconciliation differs')
    fine=pd.read_csv(out/'resolution_neighborhoods.csv.gz');fine_pairs=pd.read_csv(out/'resolution_pairs.csv');fine_summary=pd.read_csv(out/'resolution_summary.csv')
    if len(fine_summary)!=48 or fine_summary.duplicated(keys).any():raise ValueError('Resolution scenario inventory')
    if fine.duplicated(keys+['target','background']).any() or fine.background.isin(priorities).any():raise ValueError('Resolution neighborhood identity')
    constrained=fine[fine.comparison=='stage_peak_abundance_dominance']
    if (constrained.target_peak_stage!=constrained.background_peak_stage).any() or (constrained.target_peak_group!=constrained.background_peak_group).any():raise ValueError('Joint stage constraint')
    if (fine.abundance_gap>fine.abundance_caliper+1e-11).any():raise ValueError('Resolution abundance caliper')
    for key,d in fine_pairs.groupby(keys[:3]):
        sets=[set(v.target) for _,v in d.groupby('comparison')]
        if len(sets)!=3 or any(s!=sets[0] for s in sets):raise ValueError('Different target sets across resolution schemes')
    for key,d in fine_pairs.groupby(keys):
        row=fine_summary
        for k,v in zip(keys,key):row=row[row[k]==v]
        excess=d.target_logFC-d.background_mean_logFC
        observed=row.iloc[0]
        if not np.isclose(excess.mean(),observed.mean_logFC_difference,rtol=1e-10,atol=1e-11) or len(d)!=observed.matched_priorities:raise ValueError('Resolution effect reconciliation')
    return dict(scenarios=len(summary),priority_records=len(eligibility),neighborhood_members=len(edges),pairs=len(pairs),age_confirmed=source['age_p14_confirmed'])


def reproduce(out):
    with tempfile.TemporaryDirectory(prefix='stage-specificity-') as folder:
        fresh=Path(folder)
        for script in ['run_stage_specificity.py','run_stage_resolution.py','audit_stage_sources.py','plot_stage_specificity.py']:
            subprocess.run(['bash','envs/python.sh','scripts/'+script,'--output',str(fresh)],cwd=ROOT,check=True)
        check(fresh)
        for name in TABLES:
            pd.testing.assert_frame_equal(pd.read_csv(out/name),pd.read_csv(fresh/name),check_exact=False,rtol=1e-9,atol=1e-10)
        for name in JSON_FILES:
            if json.loads((out/name).read_text())!=json.loads((fresh/name).read_text()):raise ValueError('JSON reproduction differs: '+name)
        for name in FIGURES:
            if digest(out/name)!=digest(fresh/name):raise ValueError('Figure reproduction differs')
        for p in (out/'prepared').glob('*.tsv.gz'):
            pd.testing.assert_frame_equal(pd.read_csv(p,sep='\t'),pd.read_csv(fresh/'prepared'/p.name,sep='\t'))
    record=dict(status='passed',scope='Same cached processed input reconstruction and matching; not biological validation or a fresh network retrieval.',tables=TABLES,JSON=JSON_FILES,figures=FIGURES,matrices=4,rtol=1e-9,atol=1e-10)
    (out/'reproduction.json').write_text(json.dumps(record,indent=2)+'\n');return record


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--reproduce',action='store_true');a=parser.parse_args();out=ROOT/'results/stage_specificity';print(check(out))
    if a.reproduce:print(reproduce(out))
