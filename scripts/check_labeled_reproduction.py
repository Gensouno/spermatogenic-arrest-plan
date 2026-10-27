#!/usr/bin/env python3
"""Fresh-directory reproduction from three frozen original matrix windows."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def main():
    target=Path(tempfile.mkdtemp(prefix='spermatogenesis-labeled-reference-'))
    for folder in ['code','scripts','tests','envs']:
        (target/folder).mkdir()
        for p in (ROOT/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:shutil.copy2(p,target/folder/p.name)
    files=['data/labeled_reference_assets.json','data/emptydrops_label_assets.json','data/labeled_reference_contract.json','labeled_reference_snapshot/README.md']
    files += [r['path'] for r in json.loads((ROOT/'data/labeled_reference_assets.json').read_text())['assets']]
    files += [r['path'] for r in json.loads((ROOT/'data/emptydrops_label_assets.json').read_text())['assets']]
    files += list(json.loads((ROOT/'labeled_reference_snapshot/design.json').read_text())['fixed_inputs_sha256'])
    for name in set(files):
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p)
    (target/'.venv').symlink_to(ROOT/'.venv',target_is_directory=True)
    logs=ROOT/'logs/labeled_reference';logs.mkdir(parents=True,exist_ok=True)
    print('Fresh directory:',target,flush=True)
    with (logs/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_labeled_reference.py')],cwd=target,
            stdout=log,stderr=subprocess.STDOUT,check=True)
    tables=['all_group_inventory.csv','profile_evidence.csv','library_evidence.csv','cell_QC.csv.gz','type_profiles.csv.gz',
        'holdout_metrics.csv','full_fit_weights.csv','library_summary.csv','label_summary.csv','fold_summary.csv','design_conditions.csv','fold_weights.csv',
        'zygotene_overlap_metrics.csv','zygotene_overlap_weights.csv']
    for name in tables:
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/'results/labeled_reference'/name),
            pd.read_csv(target/'results/labeled_reference'/name),rtol=1e-10,atol=1e-12)
    for name in ['preparation_audit.json','design.json','emptydrops_label_audit.json','zygotene_overlap_audit.json']:
        first=json.loads((ROOT/'results/labeled_reference'/name).read_text());second=json.loads((target/'results/labeled_reference'/name).read_text())
        first.pop('elapsed_seconds',None);second.pop('elapsed_seconds',None)
        if first!=second:raise ValueError('Source/design differs')
    report={'status':'complete','tables_matching':len(tables),'source_design_JSON_matching':4,
        'scope':'Three cached original byte windows reprocessed in a fresh directory, frozen references and shared installed Python dependencies reused; no new network transfer, full-atlas reprocessing or clean installation.'}
    (ROOT/'labeled_reference_snapshot/reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print('Labelled-library reproduction passed:',len(tables),'tables and four JSON files')


if __name__=='__main__':main()
