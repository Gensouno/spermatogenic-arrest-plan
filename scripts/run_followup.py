#!/usr/bin/env python3
"""Run the isolated follow-up; network access requires --download explicitly."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    os.chdir(ROOT)
    env = os.environ.copy()
    env['R_LIBS_USER'] = str(ROOT/'.Rlib') + (os.pathsep+env['R_LIBS_USER'] if env.get('R_LIBS_USER') else '')
    env['MPLCONFIGDIR'] = str(ROOT/'.mplcache')
    steps = [[sys.executable, 'scripts/fetch_followup.py'] + (['--download'] if args.download else []),
             [sys.executable, 'code/prepare_followup.py'], ['Rscript', 'code/analyze_lrh1.R'],
             [sys.executable, 'code/summarize_followup.py'],
             [sys.executable, '-m', 'unittest', 'discover', '-s', 'tests', '-v']]
    logdir = ROOT/'logs/followup';logdir.mkdir(parents=True, exist_ok=True)
    record = {'download_requested': args.download, 'steps': []}
    for i, command in enumerate(steps):
        start = time.monotonic();print('RUN', command[1], flush=True)
        with (logdir/f'step_{i:02d}.log').open('w') as log:
            result = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        record['steps'].append({'command': command[1:], 'exit_code': result.returncode,
                                'elapsed_seconds': round(time.monotonic()-start, 3)})
        record['status'] = 'running' if result.returncode == 0 else 'failed'
        (logdir/'run.json').write_text(json.dumps(record, indent=2)+'\n')
        if result.returncode:
            raise SystemExit(f'Stopped; inspect {logdir / ("step_%02d.log" % i)}')
    record['status'] = 'complete'
    (logdir/'run.json').write_text(json.dumps(record, indent=2)+'\n')
    print('Isolated follow-up complete; previous snapshots retained.')


if __name__ == '__main__':
    main()
