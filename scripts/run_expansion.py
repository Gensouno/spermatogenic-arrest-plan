#!/usr/bin/env python3
"""Optional data expansion; preserve the frozen eight-cohort baseline."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true', help='Fetch missing metadata/assets from explicit primary sources; otherwise require cached inputs')
    args = parser.parse_args()
    review = json.loads((ROOT / 'data/expansion_review.json').read_text())
    assets = json.loads((ROOT / 'data/expansion_assets.json').read_text())['assets']
    required_baseline = ['results/all_cohort_comparison/common_logFC.csv.gz',
                         'results/stage_reference/pooled_centroid_13_reference.csv.gz']
    if any(not (ROOT / p).exists() for p in required_baseline):
        raise SystemExit('Run the baseline scripts/run_analysis.py first; baseline intermediates are required')
    if not args.download:
        missing = [r['gse'] for r in review['candidates'] if not (ROOT / f'data/raw/expansion_metadata/{r["gse"]}_family.soft.gz').exists()]
        missing += [r['path'] for r in assets if not (ROOT / r['path']).exists()]
        if missing:
            raise SystemExit('Missing cached expansion input; use --download explicitly: ' + ', '.join(missing))
    os.chdir(ROOT)
    env = os.environ.copy()
    env['R_LIBS_USER'] = str(ROOT / '.Rlib') + (os.pathsep + env['R_LIBS_USER'] if env.get('R_LIBS_USER') else '')
    env['MPLCONFIGDIR'] = str(ROOT / '.mplcache')
    steps = [[sys.executable, 'scripts/discover_expansion.py', '--accessions', *[r['gse'] for r in review['candidates']]]]
    if args.download:
        steps += [[sys.executable, 'scripts/fetch_expansion_asset.py', '--gse', r['gse'], '--url', r['url']] for r in assets]
    steps += [[sys.executable, 'code/prepare_expansion.py'], ['Rscript', 'code/de_expansion.R'],
              ['Rscript', 'code/map_sorted_expansion.R'], [sys.executable, 'code/check_sorted_expansion.py'],
              [sys.executable, 'code/summarize_expansion.py'],
              [sys.executable, '-m', 'unittest', 'discover', '-s', 'tests', '-v']]
    out = ROOT / 'logs/expansion'; out.mkdir(parents=True, exist_ok=True)
    record = {'download_requested': args.download, 'steps': [],
              'baseline_inputs_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in required_baseline}}
    for i, command in enumerate(steps):
        started = time.monotonic()
        print('RUN', command[1], flush=True)
        with (out / f'step_{i:02d}.log').open('w') as log:
            result = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        record['steps'].append({'command': command[1:], 'exit_code': result.returncode,
                                'elapsed_seconds': round(time.monotonic() - started, 3)})
        record['status'] = 'running' if result.returncode == 0 else 'failed'
        (out / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
        if result.returncode:
            raise SystemExit(f'Expansion stopped; see {out / ("step_%02d.log" % i)}')
    record['status'] = 'complete'
    (out / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
    print('Optional expansion complete; baseline outputs retained. Logs:', out)


if __name__ == '__main__':
    main()
