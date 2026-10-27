#!/usr/bin/env python3
"""Independent source diagnostic, offline by default; no reference training."""
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
    if args.download:
        subprocess.run([sys.executable,str(ROOT/'scripts/discover_leptotene.py'),'--download'],cwd=ROOT,env=env,check=True)
    subprocess.run([sys.executable,str(ROOT/'scripts/fetch_leptotene.py'),*(['--download'] if args.download else [])],cwd=ROOT,env=env,check=True)
    subprocess.run([sys.executable,str(ROOT/'code/diagnose_leptotene.py')],cwd=ROOT,env=env,check=True)


if __name__=='__main__':main()
