#!/usr/bin/env python3
"""Optional three-library author-label audit; offline unless --download supplied."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--download',action='store_true');a=p.parse_args()
    env=os.environ.copy()
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:env[key]='1'
    for script,extra in [('scripts/fetch_labeled_ranges.py',['--download'] if a.download else []),
        ('scripts/audit_emptydrops_labels.py',['--download'] if a.download else []),
        ('code/prepare_labeled_reference.py',[]),('code/check_author_zygotene.py',[]),('code/diagnose_labeled_reference.py',[])]:
        subprocess.run([sys.executable,str(ROOT/script),*extra],cwd=ROOT,env=env,check=True)


if __name__=='__main__':main()
