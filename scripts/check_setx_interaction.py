#!/usr/bin/env python3
"""Offline numerical checks and fresh-directory Setx interaction reproduction."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from run_setx_interaction import check_inputs,digest,summarize_candidates
from check_project import verify_manifest

ROOT=Path(__file__).resolve().parents[1]
TABLES=['sample_design.csv','all_gene_statistics.csv.gz','candidate_effects.csv.gz',
        'candidate_stability.csv','variant_summary.csv','runtime_R_versions.csv']
AUDITS=['design_audit.json','summary.json']
FIGURES=['fixed_priorities_interaction.png']


def check_results(root):
    check_inputs(root);out=root/'results/setx_interaction'
    stats=pd.read_csv(out/'all_gene_statistics.csv.gz')
    if stats.duplicated(['variant','gene']).any():raise ValueError('Duplicate gene test')
    summaries=pd.read_csv(out/'variant_summary.csv').set_index('variant')
    if set(summaries.index)!={'main','pseudo_0p1','pseudo_1','exclude_ambiguous_pair'}:
        raise ValueError('Incorrect sensitivity family')
    for name,d in stats.groupby('variant'):
        row=summaries.loc[name]
        if len(d)!=row.tested_genes or not np.isfinite(d[['F','P','BY','BH']].to_numpy()).all():
            raise ValueError('Incomplete or nonfinite omnibus tests')
        for method in ['by','bh']:
            if not np.allclose(false_discovery_control(d.P.to_numpy(),method=method),d[method.upper()],rtol=1e-10,atol=1e-12):
                raise ValueError('Whole-family adjustment mismatch')
            if int((d[method.upper()]<.05).sum())!=row['conditional_'+method.upper()+'05']:
                raise ValueError('Incorrect rejection summary')
        if row.design_rank!=10 or row.interaction_df!=4 or row.residual_df!=row.samples-10:
            raise ValueError('Invalid design degrees of freedom')
        if row.direct_coefficient_max_error>1e-10 or row.nested_F_max_relative_error>1e-8:
            raise ValueError('Independent coefficient/F audit failed')
    effects=pd.read_csv(out/'candidate_effects.csv.gz')
    candidates=pd.read_csv(root/'replicability_snapshot/candidate_stage_context.csv')
    stable=set(pd.read_csv(root/'replicability_snapshot/stable_descriptive_priorities.csv').gene)
    expected=summarize_candidates(effects,candidates,stable)
    pd.testing.assert_frame_equal(expected,pd.read_csv(out/'candidate_stability.csv'),rtol=1e-10,atol=1e-12)
    print('Setx design, nested-model F, full-family BY/BH and missing-aware candidate checks passed.',flush=True)


def reproduce(root):
    check_results(root)
    target=Path(tempfile.mkdtemp(prefix='setx-interaction-reproduce-'))
    names=['scripts/run_setx_interaction.py','code/analyze_setx_interaction.R',
           'data/setx_interaction_inputs.json','data/setx_interaction_contract.json',
           'envs/python.sh','envs/rscript.sh','envs/network.sh']
    names += [a['path'] for a in json.loads((root/'data/setx_interaction_inputs.json').read_text())['assets']]
    for name in names:
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/name,p)
    for name in ['.venv','.Rlib']:(target/name).symlink_to((root/name).resolve(),target_is_directory=True)
    logs=root/'logs/setx_interaction';logs.mkdir(parents=True,exist_ok=True)
    with (logs/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_setx_interaction.py')],cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    out=root/'results/setx_interaction';fresh=target/'results/setx_interaction'
    for name in TABLES:
        pd.testing.assert_frame_equal(pd.read_csv(out/name),pd.read_csv(fresh/name),rtol=1e-10,atol=1e-12)
    for name in AUDITS:
        if json.loads((out/name).read_text())!=json.loads((fresh/name).read_text()):
            raise ValueError('Reproduction JSON differs: '+name)
    for name in FIGURES:
        if digest(out/name)!=digest(fresh/name):raise ValueError('Reproduction image differs')
    # Gzip headers need not match; compare complete reconstructed numerical matrix.
    pd.testing.assert_frame_equal(pd.read_csv(out/'TPM.tsv.gz',sep='\t'),pd.read_csv(fresh/'TPM.tsv.gz',sep='\t'))
    report=dict(status='complete',tables_matching=len(TABLES),JSON_matching=len(AUDITS),figures_matching=len(FIGURES),
        original_genome_members_reprocessed=40,conditional_omnibus_fits=4,effect_only_omission_fits=40,
        scope='Fresh directory, original frozen source TAR and metadata re-read, all fits and fixed-candidate summaries regenerated. '
        'Installed dependencies and frozen gene annotation reused; no raw-read alignment, clean installation, new download, '
        'animal-independence verification or empirical biological FDR validation.')
    (out/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reproduce',action='store_true')
    args=parser.parse_args()
    if args.reproduce:reproduce(ROOT)
    else:
        check_results(ROOT)
        if (ROOT/'setx_interaction_snapshot/manifest.json').exists():
            print(verify_manifest(ROOT,'setx_interaction_snapshot'))


if __name__=='__main__':main()
