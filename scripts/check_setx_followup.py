#!/usr/bin/env python3
"""Numerical checks and offline fresh-directory reproduction for Setx follow-up."""
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
from audit_setx_sources import check_inputs

ROOT=Path(__file__).resolve().parents[1]
SOURCE_TABLES=['sample_identity.csv','ambiguous_pair.csv']
SOURCE_JSON=['summary.json','pair_attribute_audit.json','publication_evidence.json']
LENGTH_TABLES=['all_gene_statistics.csv.gz','variant_summary.csv','candidate_context.csv','candidate_summary.csv',
    'matched_universe_comparison.csv','runtime_versions.csv']


def check(root):
    check_inputs(root)
    out=root/'results/setx_length'
    summary=pd.read_csv(out/'variant_summary.csv').set_index('variant')
    for name,d in pd.read_csv(out/'all_gene_statistics.csv.gz').groupby('variant'):
        if d.gene.duplicated().any() or len(d)!=summary.loc[name,'tested_genes']:raise ValueError('Gene test denominator differs')
        if not np.isfinite(d[['F','P','BY','BH']].to_numpy()).all():raise ValueError('Nonfinite model statistics')
        for method in ['by','bh']:
            if not np.allclose(false_discovery_control(d.P.to_numpy(),method=method),d[method.upper()],rtol=1e-10,atol=1e-12):
                raise ValueError('Length-model complete-family adjustment differs')
    if summary.counts_rounded.any() or (summary.length_offset_max_error>1e-12).any():raise ValueError('Invalid estimated-count/offset handling')
    print('Source identities and full-family length-model BY/BH checks passed.',flush=True)


def reproduce(root):
    check(root);target=Path(tempfile.mkdtemp(prefix='setx-source-length-reproduction-'))
    for folder in ['scripts','code','tests','envs']:
        (target/folder).mkdir()
        for p in (root/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:shutil.copy2(p,target/folder/p.name)
    names=['data/setx_source_audit_contract.json','data/setx_source_audit_inputs.json','data/setx_length_contract.json','data/setx_length_inputs.json','data/setx_interaction_inputs.json']
    for manifest in ['data/setx_source_audit_inputs.json','data/setx_length_inputs.json','data/setx_interaction_inputs.json']:
        names += [v['path'] for v in json.loads((root/manifest).read_text())['assets']]
    names += ['docs/SETX_INTERACTION.md','data/setx_interaction_contract.json']
    for name in set(names):
        dest=target/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/name,dest)
    shutil.copytree(root/'setx_interaction_snapshot',target/'setx_interaction_snapshot',dirs_exist_ok=True)
    for n in ['.venv','.Rlib']:(target/n).symlink_to((root/n).resolve(),target_is_directory=True)
    logs=root/'logs/setx_source_audit';logs.mkdir(exist_ok=True,parents=True)
    with (logs/'fresh_run.log').open('w') as log:
        for script in ['scripts/audit_setx_sources.py','scripts/run_setx_length.py']:
            subprocess.run([sys.executable,str(target/script)],cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    for branch,tables in [('setx_source_audit',SOURCE_TABLES),('setx_length',LENGTH_TABLES)]:
        for name in tables:
            pd.testing.assert_frame_equal(pd.read_csv(root/'results'/branch/name),pd.read_csv(target/'results'/branch/name),rtol=1e-10,atol=1e-12)
    for branch,names in [('setx_source_audit',SOURCE_JSON),('setx_length',['summary.json'])]:
        for name in names:
            if json.loads((root/'results'/branch/name).read_text())!=json.loads((target/'results'/branch/name).read_text()):
                raise ValueError('Source/length JSON differs: '+branch+'/'+name)
    for name in ['expected_count.tsv.gz','effective_length.tsv.gz']:
        pd.testing.assert_frame_equal(pd.read_csv(root/'results/setx_length'/name,sep='\t'),pd.read_csv(target/'results/setx_length'/name,sep='\t'))
    report=dict(status='complete',tables_matching=len(SOURCE_TABLES)+len(LENGTH_TABLES),JSON_matching=len(SOURCE_JSON)+1,
        original_genome_members_reprocessed=40,length_offset_models_refitted=2,
        scope='Fresh-directory metadata audit and original RSEM estimate/length reconstruction; installed dependencies, frozen prior annotations and TPM comparison reused. No network rerun, raw-read alignment, clean installation or biological validation.')
    for branch in ['setx_source_audit','setx_length']:(root/'results'/branch/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reproduce',action='store_true');a=p.parse_args()
    reproduce(ROOT) if a.reproduce else check(ROOT)
