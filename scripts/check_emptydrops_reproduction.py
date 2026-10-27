#!/usr/bin/env python3
"""Fresh-directory rerun of both count branches using local frozen inputs."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

from check_emptydrops_p15 import input_rows
from export_emptydrops_p15 import TABLES,AUDITS

ROOT=Path(__file__).resolve().parents[1]


def main():
    target=Path(tempfile.mkdtemp(prefix='spermatogenesis-emptydrops-p15-'))
    for folder in ['code','scripts','tests','envs']:
        (target/folder).mkdir()
        for p in (ROOT/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:shutil.copy2(p,target/folder/p.name)
    files=['data/emptydrops_p15_assets.json','data/labeled_reference_assets.json','data/emptydrops_p15_contract.json',
        'labeled_reference_snapshot/design.json']+[r['path'] for r in input_rows(ROOT)]
    files+=list(json.loads((ROOT/'results/emptydrops_p15/design.json').read_text())['fixed_inputs_sha256'])
    for name in set(files):
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p)
    (target/'.venv').symlink_to(ROOT/'.venv',target_is_directory=True)
    logs=ROOT/'logs/emptydrops_p15';logs.mkdir(parents=True,exist_ok=True)
    print('Fresh directory:',target,flush=True)
    with (logs/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_emptydrops_p15.py')],cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    tables=TABLES+['cell_QC.csv.gz','cellranger_cell_QC.csv.gz','type_profiles.csv.gz','gene_design.csv.gz']
    for name in tables:
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/'results/emptydrops_p15'/name),
            pd.read_csv(target/'results/emptydrops_p15'/name),rtol=1e-10,atol=1e-12)
    for name in AUDITS:
        first=json.loads((ROOT/'results/emptydrops_p15'/name).read_text());second=json.loads((target/'results/emptydrops_p15'/name).read_text())
        first.pop('elapsed_seconds',None);second.pop('elapsed_seconds',None)
        if first!=second:raise ValueError('Reproduced audit/design differs')
    report={'status':'complete','tables_matching':len(tables),'source_design_JSON_matching':len(AUDITS),
        'scope':'Fresh-directory offline rerun from cached complete P15 EmptyDrops and CellRanger windows, unchanged references and shared installed dependencies; not a clean installation or full atlas reprocessing.'}
    (ROOT/'emptydrops_p15_snapshot/reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print('Fresh-directory reproduction passed:',len(tables),'tables and',len(AUDITS),'JSON files',flush=True)


if __name__=='__main__':main()
