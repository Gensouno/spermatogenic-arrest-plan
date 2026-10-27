#!/usr/bin/env python3
"""Run source reconciliation, technical sensitivity and transfer, offline by default."""
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
    stages=[('scripts/fetch_chen_umi.py',['--download'] if args.download else []),('code/prepare_chen_reference.py',[]),
        ('code/diagnose_reference_limits.py',[]),('code/check_chen_transfer.py',[]),('code/check_chen_library_holdout.py',[]),
        ('code/check_legacy_sorted_transfer.py',[]),('code/plot_reference_limits.py',[])]
    for script,extra in stages:subprocess.run([sys.executable,str(ROOT/script),*extra],cwd=ROOT,env=env,check=True)


if __name__=='__main__':main()
