#!/usr/bin/env python3
"""Fresh-directory offline refit from frozen prepared matrices and stage sources."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import pandas as pd
from check_replicability import check_inputs
from export_replicability import TABLES,AUDITS,PRIVATE_TABLES,FIGURES
from freeze_replicability_inputs import digest
ROOT=Path(__file__).resolve().parents[1]


def main():
    check_inputs(ROOT);target=Path(tempfile.mkdtemp(prefix='spermatogenesis-replicability-'))
    for folder in ['code','scripts','tests','envs']:
        (target/folder).mkdir()
        for p in (ROOT/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:shutil.copy2(p,target/folder/p.name)
    paths=[r['path'] for r in json.loads((ROOT/'data/replicability_inputs.json').read_text())['assets']]
    paths+=['data/replicability_inputs.json','data/replicability_contract.json','data/replicability_context_contract.json']
    for name in set(paths):
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p)
    for n in ['.venv','.Rlib']:(target/n).symlink_to((ROOT/n).resolve(),target_is_directory=True)
    logs=ROOT/'logs/replicability';logs.mkdir(parents=True,exist_ok=True)
    print('Fresh directory:',target,flush=True)
    with (logs/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_replicability.py')],cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    for name in TABLES+PRIVATE_TABLES:
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/'results/replicability'/name),
            pd.read_csv(target/'results/replicability'/name),rtol=1e-10,atol=1e-12)
    for name in AUDITS:
        if json.loads((ROOT/'results/replicability'/name).read_text())!=json.loads((target/'results/replicability'/name).read_text()):
            raise ValueError('Reproduced design differs: '+name)
    for name in FIGURES:
        if digest(ROOT/'results/replicability'/name)!=digest(target/'results/replicability'/name):
            raise ValueError('Reproduced figure differs: '+name)
    report={'status':'complete','tables_matching':len(TABLES+PRIVATE_TABLES),'design_JSON_matching':len(AUDITS),
        'figure_files_matching':len(FIGURES),
        'scope':'Fresh-directory offline re-fit of 8 prepared expression cohorts, all 40 primary count one-sample omissions, four-unit replicability and expression-context annotation. Six external contexts recomputed only from frozen prior fits/LOO caches, without external refits. Installed dependencies and frozen prepared input matrices reused. No raw CEL normalization, transcript annotation rebuild, clean installation, fresh download or biological validation.'}
    snap=ROOT/'replicability_snapshot';snap.mkdir(exist_ok=True)
    (snap/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
