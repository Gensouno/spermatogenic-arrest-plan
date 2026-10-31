#!/usr/bin/env python3
"""Run a separate post-hoc competitive family on four audited count contexts."""
import argparse
import json
import subprocess
import platform
import importlib.metadata
from pathlib import Path
import pandas as pd
from check_project import verify_manifest
from run_independent_contexts import ROOT, COHORTS, prepare, digest


def run(out):
    for snapshot in ['replicability_snapshot', 'independent_contexts_snapshot']:
        verify_manifest(ROOT, snapshot)
    out.mkdir(parents=True, exist_ok=True)
    prepared=out/'prepared';prepared.mkdir(exist_ok=True)
    prepare(ROOT, prepared)
    subprocess.run(['bash','envs/rscript.sh','code/analyze_competitive_set.R',str(out)],cwd=ROOT,check=True)
    results=pd.read_csv(out/'camera_results.csv')
    summary=dict(scope='Post-hoc two-sided competitive association; not primary replication or biological validation',
        family_tests=len(results),correlation_policies=['estimated','fixed_0p01'],
        fixed_priorities=36,primary_discoveries_changed=False,
        results={p:dict(contexts=int(len(d)),direction_down=int(d.Direction.eq('Down').sum()),
            BY_below_0p05=d.loc[d.BY.lt(.05),'cohort'].tolist(),BH_below_0p05=d.loc[d.BH.lt(.05),'cohort'].tolist())
            for p,d in results.groupby('policy',sort=False)},
        contracts={'data/competitive_set_contract.json':digest(ROOT/'data/competitive_set_contract.json')},
        inherited_input_specification={'data/independent_contexts_inputs.json':digest(ROOT/'data/independent_contexts_inputs.json')},
        random_seed=None)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    runtime={'Python':platform.python_version(),**{p:importlib.metadata.version(p) for p in ['numpy','pandas','matplotlib']}}
    (out/'runtime_Python_versions.json').write_text(json.dumps(runtime,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results/competitive_set')
    run(parser.parse_args().output.resolve())
