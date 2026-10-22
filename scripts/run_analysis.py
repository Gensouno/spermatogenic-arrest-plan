#!/usr/bin/env python3
"""Run audited analysis stages in dependency order, with logs and source hashes."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
STAGES = [
    ('counts', 'R', 'code/prepare_counts_and_de.R'),
    ('legacy_annotation', 'python', 'code/prepare_legacy_annotation.py'),
    ('extended_import', 'python', 'code/import_extended_rnaseq.py'),
    ('extended_de', 'R', 'code/de_extended_rnaseq.R'),
    ('arrays_import', 'python', 'code/import_microarrays.py'),
    ('arrays_de', 'R', 'code/de_microarrays.R'),
    ('sample_metadata', 'python', 'code/verify_count_metadata.py'),
    ('reference_aggregate', 'python', 'code/build_reference_batches.py'),
    ('reference_calibrate', 'python', 'code/calibrate_reference_batches.py'),
    ('mixture_audit', 'python', 'code/audit_mixture.py'),
    ('atlas_metadata', 'python', 'code/audit_reference_metadata.py'),
    ('effect_compare', 'python', 'code/compare_rnaseq_effects.py', '--arrays'),
    ('composition_explore', 'python', 'code/explore_count_cohorts.py', '--extended-counts'),
    ('composition_crossfit', 'python', 'code/crossfit_composition_prediction.py'),
    ('robustness', 'R', 'code/check_effect_robustness.R'),
    ('robustness_summary', 'python', 'code/summarize_robustness.py'),
    ('external_sorted', 'python', 'code/check_external_sorted_reference.py'),
    ('fine_stage_reference', 'python', 'code/check_stage_reference.py'),
    ('tests', 'python', '-m', 'unittest', 'discover', '-s', 'tests', '-v'),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    names = [stage[0] for stage in STAGES]
    parser.add_argument('--from-stage', choices=names, default=names[0])
    parser.add_argument('--to-stage', choices=names, default=names[-1])
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--skip-asset-check', action='store_true', help='For debugging existing outputs only; recorded in manifest')
    args = parser.parse_args()
    if names.index(args.from_stage) > names.index(args.to_stage):
        parser.error('from-stage must not follow to-stage')
    selected = STAGES[names.index(args.from_stage):names.index(args.to_stage) + 1]
    if args.list:
        for name, runtime, *command in selected:
            print(name, runtime, *command)
        return
    os.chdir(ROOT)
    for directory in ['data/bulk', 'results/de', 'results/import_qc', 'logs/analysis']:
        (ROOT / directory).mkdir(parents=True, exist_ok=True)
    rscript = shutil.which('Rscript')
    if any(s[1] == 'R' for s in selected) and not rscript:
        raise RuntimeError('Rscript not found in PATH')
    env = os.environ.copy()
    env['R_LIBS_USER'] = str(ROOT / '.Rlib') + (os.pathsep + env['R_LIBS_USER'] if env.get('R_LIBS_USER') else '')
    env['MPLCONFIGDIR'] = str(ROOT / '.mplcache')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    out = ROOT / 'logs/analysis' / stamp
    out.mkdir()
    record = {'started_utc': stamp, 'python': sys.executable, 'python_version': sys.version,
        'Rscript': rscript, 'asset_check_skipped': args.skip_asset_check, 'stages': [],
        'curation_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [ROOT / 'data/samples.tsv', ROOT / 'data/analysis_assets.json']},
        'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for directory in ['code', 'scripts', 'tests'] for p in (ROOT / directory).glob('*')
            if p.is_file() and p.suffix in ['.py', '.R']}}
    status_path = out / 'run.json'
    def save():
        status_path.write_text(json.dumps(record, indent=2) + '\n')
    save()
    with (out / 'python_packages.json').open('w') as package_file:
        from importlib.metadata import distributions
        json.dump(sorted([(d.metadata['Name'], d.version) for d in distributions()]), package_file, indent=2)
    if not args.skip_asset_check:
        subprocess.run([sys.executable, str(ROOT / 'scripts/data_assets.py'), '--verify'], check=True)
    for name, runtime, *command in selected:
        executable = sys.executable if runtime == 'python' else rscript
        row = {'name': name, 'command': [executable, *command], 'status': 'running'}
        record['stages'].append(row)
        save()
        print(f'RUN {name}', flush=True)
        started = time.monotonic()
        with (out / f'{name}.log').open('w') as log:
            process = subprocess.Popen([executable, *command], cwd=ROOT, env=env,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                print(line, end='', flush=True)
                log.write(line)
                log.flush()
            code = process.wait()
        row.update(status='complete' if code == 0 else 'failed', exit_code=code,
                   elapsed_seconds=round(time.monotonic() - started, 3))
        save()
        if code:
            raise SystemExit(f'Stopped at {name}; see {out / (name + ".log")}')
    record['status'] = 'complete'
    save()
    print('Completed selected stages; provenance:', status_path)


if __name__ == '__main__':
    main()
