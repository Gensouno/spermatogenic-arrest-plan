#!/usr/bin/env python3
"""Fresh-directory offline reprocessing of the two new RNA study sources."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

from export_leptotene import TABLES,AUDITS

ROOT=Path(__file__).resolve().parents[1]


def main():
    target=Path(tempfile.mkdtemp(prefix='spermatogenesis-leptotene-'))
    for folder in ['code','scripts','tests','envs']:
        (target/folder).mkdir()
        for p in (ROOT/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:shutil.copy2(p,target/folder/p.name)
    files=['data/leptotene_assets.json','data/leptotene_paper_asset.json','data/leptotene_contract.json',
        'data/labeled_reference_assets.json','labeled_reference_snapshot/design.json']
    files += [r['path'] for r in json.loads((ROOT/'data/leptotene_assets.json').read_text())['assets']]
    files += list(json.loads((ROOT/'results/leptotene_followup/design.json').read_text())['fixed_inputs_sha256'])
    for name in set(files):
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p)
    (target/'.venv').symlink_to(ROOT/'.venv',target_is_directory=True)
    logs=ROOT/'logs/leptotene_followup';logs.mkdir(parents=True,exist_ok=True)
    print('Fresh directory:',target,flush=True)
    with (logs/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_leptotene.py')],cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    tables=TABLES+['common_gene_counts.csv.gz']
    for name in tables:
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/'results/leptotene_followup'/name),pd.read_csv(target/'results/leptotene_followup'/name),rtol=1e-10,atol=1e-12)
    for name in AUDITS:
        if json.loads((ROOT/'results/leptotene_followup'/name).read_text())!=json.loads((target/'results/leptotene_followup'/name).read_text()):raise ValueError('Source/design JSON differs')
    report={'status':'complete','tables_matching':len(tables),'source_design_JSON_matching':len(AUDITS),
        'scope':'Fresh-directory offline recomputation of 10 source RNA profiles from frozen raw counts/metadata, original paper XML, fixed gene mapping and unchanged references; shared installed dependencies, not a clean installation or new network-transfer validation.'}
    (ROOT/'leptotene_snapshot/reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print('Reproduced',len(tables),'tables and',len(AUDITS),'JSON files',flush=True)


if __name__=='__main__':main()
