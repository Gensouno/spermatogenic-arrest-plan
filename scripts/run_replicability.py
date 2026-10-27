#!/usr/bin/env python3
"""Offline source refit, full count LOO and independent context analysis."""
from pathlib import Path
import os
import subprocess
import sys
from check_replicability import check_inputs
ROOT=Path(__file__).resolve().parents[1]


def main():
    print('Verified inputs:',check_inputs(ROOT),flush=True)
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    for source in ['audit_replicability_statistics.R','audit_replicability_LOO.R','analyze_replicability.py',
                   'audit_replicability_adjustments.R','annotate_replicability.py','check_replicability_context.py','plot_replicability.py']:
        command=['bash','envs/rscript.sh','code/'+source] if source.endswith('.R') else [sys.executable,'code/'+source]
        print('Running',source,flush=True);subprocess.run(command,cwd=ROOT,env=env,check=True)


if __name__=='__main__':main()
