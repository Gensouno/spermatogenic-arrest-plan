#!/usr/bin/env python3
"""Verify full model families and reproduce Setx count omissions offline."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from run_setx_stability import check_inputs,REQUIRED

ROOT=Path(__file__).resolve().parents[1]
TABLES=['candidate_omissions.csv.gz','candidate_stability.csv','omission_models.csv','previous_fit_agreement.csv',
    'matched_diagnostics.csv.gz','matched_count_normalization.csv','matched_design_audit.csv','matched_summary.csv','runtime_versions.csv']
JSON_FILES=['summary.json','matched_threshold_transitions.json']
PRIVATE_TABLES=['all_gene_omissions.csv.gz']
FIGURES=['stability_and_matched_models.png']


def check(root):
    check_inputs(root);out=root/'results/setx_stability'
    models=pd.read_csv(out/'omission_models.csv').set_index('variant')
    if set(models.index)!=set(REQUIRED) or models.index.duplicated().any():raise ValueError('Incomplete omission design')
    if models.counts_rounded.any() or not (models.design_rank==10).all() or not (models.interaction_df==4).all():
        raise ValueError('Count/contrast handling differs')
    d=pd.read_csv(out/'all_gene_omissions.csv.gz')
    if set(d.variant)!=set(REQUIRED) or d.duplicated(['variant','gene']).any():raise ValueError('Incomplete original gene fits')
    for name,part in d.groupby('variant'):
        row=models.loc[name]
        if len(part)!=row.tested_genes:raise ValueError('Omission test family differs')
        if name.startswith('omit_') and row.samples!=39:raise ValueError('Incorrect omission sample count')
        verify_adjustment(part)
        if (part.BY<.05).sum()!=row.conditional_BY05 or (part.BH<.05).sum()!=row.conditional_BH05:raise ValueError('Incorrect family summary')
    diagnostics=pd.read_csv(out/'matched_diagnostics.csv.gz')
    names={'matched_TPM_limma','matched_count_length_offset','matched_count_zero_offset_counterfactual'}
    if set(diagnostics.variant)!=names or diagnostics.duplicated(['variant','gene']).any():raise ValueError('Incorrect matched models')
    sets=[]
    for name,part in diagnostics.groupby('variant'):
        sets.append(set(part.gene));verify_adjustment(part)
    if not all(s==sets[0] for s in sets):raise ValueError('Matched test universes differ')
    audit=pd.read_csv(out/'matched_design_audit.csv').iloc[0]
    if audit.common_tested_genes!=len(sets[0]) or audit.count_library_normalization_max_error>1e-12:
        raise ValueError('Matched normalization/denominator audit failed')
    if (pd.read_csv(out/'previous_fit_agreement.csv').maximum_absolute_error>1e-9).any():raise ValueError('Previous frozen fits differ')
    print('42 omission families, 3 matched families, normalization and frozen-fit checks passed.',flush=True)


def verify_adjustment(d):
    if not np.isfinite(d[['F','P','BY','BH','delta_P18']].to_numpy()).all():raise ValueError('Nonfinite statistics')
    for method in ['by','bh']:
        if not np.allclose(false_discovery_control(d.P.to_numpy(),method=method),d[method.upper()],rtol=1e-10,atol=1e-12):
            raise ValueError('Full-family '+method+' differs')


def reproduce(root):
    check(root);target=Path(tempfile.mkdtemp(prefix='setx-stability-reproduction-'))
    # Include sources/specifications bound by the two reused frozen branches.
    paths={'data/setx_stability_inputs.json','data/setx_stability_contract.json'}
    for folder in ['code','scripts','tests','envs']:
        for p in (root/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:paths.add(str(p.relative_to(root)))
    paths.update(a['path'] for a in json.loads((root/'data/setx_stability_inputs.json').read_text())['assets'])
    for branch in ['setx_interaction_snapshot','setx_length_snapshot']:
        manifest=json.loads((root/branch/'manifest.json').read_text())
        paths.add(branch+'/manifest.json');paths.update(branch+'/'+a['path'] for a in manifest['files'])
        for field in ['source_sha256','input_specification_sha256']:paths.update(manifest[field])
    for name in paths:
        dest=target/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/name,dest)
    for name in ['.venv','.Rlib']:(target/name).symlink_to((root/name).resolve(),target_is_directory=True)
    logs=root/'logs/setx_stability';logs.mkdir(exist_ok=True,parents=True)
    with (logs/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_setx_stability.py')],cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    original=root/'results/setx_stability';fresh=target/'results/setx_stability'
    for name in TABLES+PRIVATE_TABLES:
        pd.testing.assert_frame_equal(pd.read_csv(original/name),pd.read_csv(fresh/name),rtol=1e-9,atol=1e-11)
    for name in JSON_FILES:
        if json.loads((original/name).read_text())!=json.loads((fresh/name).read_text()):raise ValueError('JSON differs: '+name)
    for name in FIGURES:
        if hashlib.sha256((original/name).read_bytes()).hexdigest()!=hashlib.sha256((fresh/name).read_bytes()).hexdigest():
            raise ValueError('Figure differs: '+name)
    for name in ['TPM','expected_count','effective_length']:
        pd.testing.assert_frame_equal(pd.read_csv(original/(name+'.tsv.gz'),sep='\t'),pd.read_csv(fresh/(name+'.tsv.gz'),sep='\t'))
    report=dict(status='complete',tables_matching=len(TABLES+PRIVATE_TABLES),JSON_matching=len(JSON_FILES),
        figures_matching=len(FIGURES),
        original_genome_members_reprocessed=40,count_omission_suite_fits=42,matched_universe_fits=3,
        scope='Original processed RSEM TAR re-read and all fits refit in a fresh directory with installed dependencies and frozen prior annotation/TPM comparison. No network rerun, clean installation, raw-read alignment or biological validation.')
    (original/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reproduce',action='store_true');args=p.parse_args()
    reproduce(ROOT) if args.reproduce else check(ROOT)
