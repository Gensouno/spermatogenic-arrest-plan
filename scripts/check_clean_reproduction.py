#!/usr/bin/env python3
"""Re-run in a new directory using only frozen inputs and shared dependencies.

This tests hidden file dependencies, NOT fresh installation or validation of
scientific conclusions. The temporary directory is retained for inspection.
"""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    dest = Path(tempfile.mkdtemp(prefix='spermatogenesis-clean-'))
    out = ROOT / 'results/reproducibility'
    out.mkdir(parents=True, exist_ok=True)
    record = {'clean_directory': str(dest), 'scope': 'Fresh outputs and intermediates, same installed Python/R dependency libraries; NOT a fresh dependency installation', 'status': 'running'}
    record_path = out / 'clean_run.json'
    record_path.write_text(json.dumps(record, indent=2) + '\n')
    print('Clean directory:', dest, flush=True)
    for directory in ['code', 'scripts', 'tests', 'envs']:
        (dest / directory).mkdir()
        for source in (ROOT / directory).iterdir():
            if source.is_file() and source.suffix in ['.py', '.R', '.sh', '.txt']:
                shutil.copy2(source, dest / directory / source.name)
    assets = json.loads((ROOT / 'data/analysis_assets.json').read_text())['assets']
    for relative in ['data/samples.tsv', 'data/analysis_assets.json'] + [a['path'] for a in assets]:
        target = dest / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, target)
    for library in ['.venv', '.Rlib']:
        (dest / library).symlink_to(ROOT / library, target_is_directory=True)
    start = time.monotonic()
    with (out / 'clean_run.log').open('w') as log:
        process = subprocess.Popen([sys.executable, str(dest / 'scripts/run_analysis.py')],
            cwd=dest, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            log.write(line)
            log.flush()
            if line.startswith(('RUN ', 'Completed ')):
                print(line, end='', flush=True)
        code = process.wait()
    record.update(exit_code=code, elapsed_seconds=round(time.monotonic() - start, 3))
    if code:
        record['status'] = 'failed'
        record_path.write_text(json.dumps(record, indent=2) + '\n')
        raise SystemExit(f'Clean pipeline failed; see {out / "clean_run.log"}')
    tables = [f'results/de/{gse}_de.csv' for gse in ['GSE236388', 'GSE283904', 'GSE68283', 'GSE68284', 'GSE101511', 'GSE294209', 'GSE38096', 'GSE42601']]
    tables += ['results/all_cohort_comparison/common_gene_spearman.csv',
        'results/robustness/sensitivity_metrics.csv',
        'results/robustness/pairwise_sensitivity.csv',
        'results/composition_prediction/prediction_metrics.csv',
        'results/composition_prediction/covariance_terms.csv',
        'results/external_sorted_reference/fitted_RNA_weights.csv',
        'results/stage_reference/fitted_weights.csv',
        'results/stage_reference/crossfit_metrics.csv',
        'results/mixture_audit/calibration_summary.csv',
        'results/mixture_audit/real_heldout.csv']
    comparisons = []
    try:
        for relative in tables:
            expected, actual = pd.read_csv(ROOT / relative), pd.read_csv(dest / relative)
            pd.testing.assert_frame_equal(actual, expected, rtol=1e-10, atol=1e-12)
            comparisons.append({'path': relative, 'rows': len(actual), 'matches': True})
    except Exception:
        record.update(status='comparison_failed', compared_tables=comparisons)
        record_path.write_text(json.dumps(record, indent=2) + '\n')
        raise
    record.update(status='complete', compared_tables=comparisons,
                  input_files_copied=len(assets), directory_retained=True)
    record_path.write_text(json.dumps(record, indent=2) + '\n')
    print('Clean run and all', len(comparisons), 'table comparisons passed:', record_path)


if __name__ == '__main__':
    main()
