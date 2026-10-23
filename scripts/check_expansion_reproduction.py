#!/usr/bin/env python3
"""Reproduce only the optional expansion from frozen baseline artifacts and raw inputs."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    target = Path(tempfile.mkdtemp(prefix='spermatogenesis-expansion-'))
    review = json.loads((ROOT / 'data/expansion_review.json').read_text())
    assets = json.loads((ROOT / 'data/expansion_assets.json').read_text())['assets']
    for folder in ['code', 'scripts', 'tests', 'envs']:
        (target / folder).mkdir()
        for p in (ROOT / folder).iterdir():
            if p.is_file() and p.suffix in ['.py', '.R', '.sh', '.txt']:
                shutil.copy2(p, target / folder / p.name)
    files = ['data/expansion_review.json', 'data/expansion_assets.json',
             'results/all_cohort_comparison/common_logFC.csv.gz',
             'data/raw/ref/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx']
    files += [r['path'] for r in assets]
    files += [f'data/raw/expansion_metadata/{r["gse"]}_family.soft.gz' for r in review['candidates']]
    files += [f'results/stage_reference/pooled_{name}_reference.csv.gz' for name in ['original_4', 'centroid_4', 'centroid_13']]
    for relative in files:
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, destination)
    for library in ['.venv', '.Rlib']:
        (target / library).symlink_to(ROOT / library, target_is_directory=True)
    print('Expansion clean directory:', target, flush=True)
    logdir = ROOT / 'logs/expansion'; logdir.mkdir(parents=True, exist_ok=True)
    with (logdir / 'fresh_run.log').open('w') as log:
        subprocess.run([sys.executable, str(target / 'scripts/run_expansion.py')],
                       cwd=target, stdout=log, stderr=subprocess.STDOUT, check=True)
    tables = [f'de/{c}_de.csv' for c in ['GSE259242_P9', 'GSE117706_P2', 'GSE117706_P7', 'GSE235777_P45']]
    tables += ['cohort_summary.csv', 'sample_sensitivity.csv', 'fixed_core_gene_correlations.csv',
               'cross_cohort_leave_one_out.csv', 'GSE232415_file_identity_audit.csv',
               'unresolved_count_units.csv', 'sorted_validation/fitted_weights.csv', 'sorted_validation/summary.csv']
    for relative in tables:
        before = pd.read_csv(ROOT / 'results/data_expansion' / relative)
        after = pd.read_csv(target / 'results/data_expansion' / relative)
        pd.testing.assert_frame_equal(before, after, rtol=1e-10, atol=1e-12)
    report = {'status': 'complete', 'tables_matching': len(tables), 'compared_tables': tables,
              'scope': 'New expansion outputs from copied frozen baseline effect/reference artifacts and raw expansion inputs; shared installed dependencies, no network or fresh installation.',
              'temporary_directory': str(target), 'retained': True}
    (logdir / 'fresh_run.json').write_text(json.dumps(report, indent=2) + '\n')
    public = {k: v for k, v in report.items() if k != 'temporary_directory'}
    (ROOT / 'expansion_snapshot/reproduction.json').write_text(json.dumps(public, indent=2) + '\n')
    print('Expansion reproduction passed:', len(tables), 'tables; local record:', logdir / 'fresh_run.json')


if __name__ == '__main__':
    main()
