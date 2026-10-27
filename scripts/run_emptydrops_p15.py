#!/usr/bin/env python3
"""Complete-P15 diagnostic, offline by default; fixed count sources only."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--download',action='store_true');args=parser.parse_args()
    env=os.environ.copy()
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:env[key]='1'
    for script,extra in [('scripts/fetch_emptydrops_p15.py',['--download'] if args.download else []),
        ('scripts/check_emptydrops_p15.py',['--inputs-only']),('code/prepare_emptydrops_p15.py',[]),('code/diagnose_emptydrops_p15.py',[])]:
        subprocess.run([sys.executable,str(ROOT/script),*extra],cwd=ROOT,env=env,check=True)


if __name__=='__main__':main()
