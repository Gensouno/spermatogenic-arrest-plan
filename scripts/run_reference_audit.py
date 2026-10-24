#!/usr/bin/env python3
"""Run bounded-memory age and fixed-reference diagnostics, offline by default."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--download',action='store_true');args=p.parse_args()
    env=os.environ.copy()
    for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:env[name]='1'
    env['MPLCONFIGDIR']=str(ROOT/'.mplcache')
    steps=[[sys.executable,'scripts/fetch_age_reference.py']+(['--download'] if args.download else []),
        [sys.executable,'code/prepare_age_reference.py'],[sys.executable,'code/diagnose_reference_transfer.py'],
        [sys.executable,'-m','unittest','discover','-s','tests','-v']]
    out=ROOT/'logs/reference_audit';out.mkdir(parents=True,exist_ok=True)
    record={'download_requested':args.download,'BLAS_threads':1,'steps':[]}
    for i,cmd in enumerate(steps):
        start=time.monotonic();print('RUN',cmd[1],flush=True)
        with (out/f'step_{i:02d}.log').open('w') as log:
            result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        record['steps'].append({'command':cmd[1:],'exit_code':result.returncode,'seconds':round(time.monotonic()-start,2)})
        record['status']='running' if result.returncode==0 else 'failed'
        (out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
        if result.returncode:raise SystemExit('Audit stopped; inspect '+str(out/f'step_{i:02d}.log'))
    record['status']='complete';(out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Reference audit complete; prior analyses retained.')


if __name__=='__main__':main()
