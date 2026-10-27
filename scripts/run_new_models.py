#!/usr/bin/env python3
"""Offline preparation, exploratory fits and fixed-list context summaries."""
from pathlib import Path
import os
import subprocess
import sys
from check_new_models import check_inputs
ROOT=Path(__file__).resolve().parents[1]


def main():
    print('Frozen inputs:',check_inputs(ROOT),flush=True)
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    for n in ['prepare_new_models.py','analyze_new_models.R','summarize_new_models.py']:
        command=['bash','envs/rscript.sh','code/'+n] if n.endswith('.R') else [sys.executable,'code/'+n]
        print('Running',n,flush=True);subprocess.run(command,cwd=ROOT,env=env,check=True)


if __name__=='__main__':main()
