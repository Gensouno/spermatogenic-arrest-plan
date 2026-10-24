#!/usr/bin/env python3
"""Fresh-directory follow-up reproduction with shared installed dependencies."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    target = Path(tempfile.mkdtemp(prefix='spermatogenesis-followup-'))
    for folder in ['code','scripts','tests','envs']:
        (target/folder).mkdir()
        for p in (ROOT/folder).iterdir():
            if p.is_file() and p.suffix in ['.py','.R','.sh','.txt']:
                shutil.copy2(p, target/folder/p.name)
    assets = json.loads((ROOT/'data/followup_assets.json').read_text())['assets']
    files = ['data/followup_assets.json','data/followup_contract.json',
        'results/all_cohort_comparison/common_logFC.csv.gz',
        'data/raw/ref/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx']
    files += [r['path'] for r in assets]
    for name in files:
        p = target/name;p.parent.mkdir(parents=True, exist_ok=True);shutil.copy2(ROOT/name,p)
    for name in ['.venv','.Rlib']:
        (target/name).symlink_to(ROOT/name, target_is_directory=True)
    logdir=ROOT/'logs/followup';logdir.mkdir(parents=True,exist_ok=True)
    print('Fresh directory:',target,flush=True)
    with (logdir/'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'scripts/run_followup.py')],cwd=target,
                       stdout=log,stderr=subprocess.STDOUT,check=True)
    tables = ['cohort_summary.csv','Lrh1_samples.csv','Ppp2ca_unit_audit.csv',
        'batch_sensitivity.csv','method_sensitivity.csv','leave_one_out.csv','PCA.csv','baseline_correlations.csv',
        'baseline_ranges.csv','paper_direction_checks.csv']
    tables += [f'de/{c}.csv' for c in ['Lrh1_P1_NanosCre','Lrh1_P3_DhhCre',
                                     'Lrh1_P1_NanosCre_unadjusted','Lrh1_P1_NanosCre_B2_only',
                                     'Lrh1_P1_NanosCre_edgeRQL','Lrh1_P3_DhhCre_edgeRQL']]
    for name in tables:
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/'results/followup'/name),
            pd.read_csv(target/'results/followup'/name),rtol=1e-10,atol=1e-12)
    report={'status':'complete','tables_matching':len(tables),'compared_tables':tables,
        'scope':'Cached raw follow-up inputs and frozen baseline effect matrix; shared installed dependencies, no download or fresh installation.'}
    (ROOT/'followup_snapshot/reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
    (logdir/'fresh_run.json').write_text(json.dumps({**report,'retained_directory':str(target)},indent=2)+'\n')
    print('Follow-up reproduction passed:',len(tables),'tables')


if __name__ == '__main__':
    main()
