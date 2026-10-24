#!/usr/bin/env python3
"""Fresh-directory reference audit with frozen prepared H3mm15/reference inputs."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def main():
    target=Path(tempfile.mkdtemp(prefix='spermatogenesis-reference-audit-'))
    for folder in ['code','scripts','tests','envs']:
        (target/folder).mkdir()
        for p in (ROOT/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:shutil.copy2(p,target/folder/p.name)
    files=['data/age_reference_assets.json','data/reference_audit_contract.json',
        'data/raw/ref/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx']
    files += [r['path'] for r in json.loads((ROOT/'data/age_reference_assets.json').read_text())['assets']]
    files += list(json.loads((ROOT/'reference_audit_snapshot/design.json').read_text())['fixed_inputs_sha256'])
    for name in set(files):
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p)
    (target/'.venv').symlink_to(ROOT/'.venv',target_is_directory=True)
    logs=ROOT/'logs/reference_audit';logs.mkdir(parents=True,exist_ok=True)
    print('Fresh directory:',target,flush=True)
    with (logs/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_reference_audit.py')],cwd=target,
                       stdout=log,stderr=subprocess.STDOUT,check=True)
    tables=['age_reference/'+n for n in ['author_sample_table.csv','library_QC.csv','profile_evidence.csv','library_profiles.csv.gz']]
    tables += ['reference_transfer/'+n for n in ['holdout_metrics.csv','fold_weights.csv','design_conditions.csv',
              'sorted_class_matches.csv','weight_ranges.csv','fine_vs_broad.csv','gene_folds.csv']]
    for name in tables:
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/'results'/name),pd.read_csv(target/'results'/name),rtol=1e-10,atol=1e-12)
    for name in ['age_reference/source_audit.json','reference_transfer/design.json']:
        a=json.loads((ROOT/'results'/name).read_text());b=json.loads((target/'results'/name).read_text())
        a.pop('elapsed_seconds',None);b.pop('elapsed_seconds',None)
        if a!=b:raise ValueError('Source/design JSON changed: '+name)
    report={'status':'complete','tables_matching':len(tables),'design_JSON_matching':2,'compared_tables':tables,
        'scope':'Raw GSE121904 matrix/source files reprocessed; frozen previously prepared H3mm15 matrix and reference tables reused; shared Python dependencies; no network or clean installation.'}
    (ROOT/'reference_audit_snapshot/reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logs/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print('Reference audit reproduction passed:',len(tables),'tables and two design/source JSON files')


if __name__=='__main__':main()
