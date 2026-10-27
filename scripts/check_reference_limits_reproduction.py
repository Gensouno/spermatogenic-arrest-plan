#!/usr/bin/env python3
"""Fresh-directory offline rerun from original Chen matrices and frozen prior inputs."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

from export_reference_limits import TABLES,AUDITS
from check_reference_limits import check_inputs

ROOT=Path(__file__).resolve().parents[1]


def main():
    check_inputs(ROOT);target=Path(tempfile.mkdtemp(prefix='spermatogenesis-reference-limits-'))
    for folder in ['code','scripts','tests','envs']:
        (target/folder).mkdir()
        for p in (ROOT/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:shutil.copy2(p,target/folder/p.name)
    paths=[r['path'] for r in json.loads((ROOT/'data/reference_limits_inputs.json').read_text())['assets']]
    paths+=[r['path'] for r in json.loads((ROOT/'data/reference_limits_sorted_inputs.json').read_text())['assets']]
    paths+=['data/reference_limits_inputs.json','data/reference_limits_sorted_inputs.json','data/reference_limits_contract.json','data/chen_umi_assets.json','data/leptotene_assets.json']
    for name in set(paths):
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p)
    (target/'.venv').symlink_to(ROOT/'.venv',target_is_directory=True)
    logs=ROOT/'logs/reference_limits';logs.mkdir(parents=True,exist_ok=True)
    print('Fresh directory:',target,flush=True)
    with (logs/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_reference_limits.py')],cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    tables=TABLES+['chen_stage_profiles.csv.gz','chen_batch_prototypes.csv.gz']
    for name in tables:
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/'results/reference_limits'/name),pd.read_csv(target/'results/reference_limits'/name),rtol=1e-10,atol=1e-12)
    for name in AUDITS:
        if json.loads((ROOT/'results/reference_limits'/name).read_text())!=json.loads((target/'results/reference_limits'/name).read_text()):raise ValueError('Reproduced JSON differs: '+name)
    report={'status':'complete','tables_matching':len(tables),'source_design_JSON_matching':len(AUDITS),
        'scope':'Fresh-directory offline reprocessing of 34 original Chen UMI matrices and author workbooks; technical sensitivity, candidate transfer and source-library holdout recomputed. Frozen earlier prepared P15 profiles and reference class sums, prior sorted counts and installed dependencies reused. Not a full main-project rerun, clean installation or fresh network test.'}
    snap=ROOT/'reference_limits_snapshot';snap.mkdir(exist_ok=True)
    (snap/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print('Fresh-directory reproduction matched',len(tables),'tables and',len(AUDITS),'JSON files',flush=True)


if __name__=='__main__':main()
