#!/usr/bin/env python3
"""Fresh-directory offline reprocessing from original two downloaded assets."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import pandas as pd
from check_new_models import check_inputs
from export_new_models import TABLES,PRIVATE_TABLES,AUDITS,FIGURES
from freeze_replicability_inputs import digest
ROOT=Path(__file__).resolve().parents[1]


def main():
    check_inputs(ROOT);target=Path(tempfile.mkdtemp(prefix='spermatogenesis-new-models-'))
    for folder in ['code','scripts','tests','envs']:
        (target/folder).mkdir()
        for p in (ROOT/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:shutil.copy2(p,target/folder/p.name)
    paths=[r['path'] for r in json.loads((ROOT/'data/new_model_inputs.json').read_text())['assets']]
    paths+=['data/new_model_inputs.json','data/new_model_audit_contract.json']
    for name in set(paths):
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p)
    for n in ['.venv','.Rlib']:(target/n).symlink_to((ROOT/n).resolve(),target_is_directory=True)
    logs=ROOT/'logs/new_models';logs.mkdir(parents=True,exist_ok=True)
    print('Fresh directory:',target,flush=True)
    with (logs/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_new_models.py')],cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    for n in TABLES+PRIVATE_TABLES:
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/'results/new_models'/n),pd.read_csv(target/'results/new_models'/n),rtol=1e-10,atol=1e-12)
    for n in AUDITS:
        if json.loads((ROOT/'results/new_models'/n).read_text())!=json.loads((target/'results/new_models'/n).read_text()):raise ValueError('New source design differs: '+n)
    for n in FIGURES:
        if digest(ROOT/'results/new_models'/n)!=digest(target/'results/new_models'/n):raise ValueError('New source figure differs: '+n)
    report={'status':'complete','tables_matching':len(TABLES+PRIVATE_TABLES),'design_JSON_matching':len(AUDITS),'figure_files_matching':len(FIGURES),
        'scope':'Original Eef1g processed gene-count file and Setx TAR genome profiles re-read, all 6 + 40 omissions re-fit, 10 TPM pseudocount alternatives and fixed-list contexts regenerated. Frozen prior candidate tables and installed dependencies reused. No raw-read alignment, fresh network test, clean installation, source DESeq2 reproduction or animal-independence verification.'}
    snap=ROOT/'new_models_snapshot';snap.mkdir(exist_ok=True)
    (snap/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
