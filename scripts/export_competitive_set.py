#!/usr/bin/env python3
"""Freeze the separate competitive branch only after successful reproduction."""
import json
import shutil
from check_competitive_set import TABLES,JSON_FILES,TEXT_FILES,FIGURES,check
from run_competitive_set import ROOT,digest
from export_reference_limits import local_source_closure


def main():
    out=ROOT/'results/competitive_set';check(out)
    if json.loads((out/'reproduction.json').read_text())['status']!='passed':raise ValueError('Reproduction required')
    snap=ROOT/'competitive_set_snapshot'
    if (snap/'manifest.json').exists():raise ValueError('Snapshot frozen; use a new branch name')
    snap.mkdir(exist_ok=True)
    for name in TABLES+JSON_FILES+TEXT_FILES+FIGURES+['reproduction.json']:shutil.copy2(out/name,snap/name)
    sources=local_source_closure(ROOT,['scripts/run_competitive_set.py','scripts/check_competitive_set.py',
        'scripts/plot_competitive_set.py','scripts/export_competitive_set.py','tests/test_competitive_set.py'])
    sources+=['code/analyze_competitive_set.R']
    specs=['data/competitive_set_contract.json','data/independent_contexts_inputs.json',
        'docs/COMPETITIVE_SET.md','envs/python.sh','envs/rscript.sh','envs/network.sh',
        'replicability_snapshot/manifest.json','independent_contexts_snapshot/manifest.json']
    manifest=dict(scope='Post-hoc competitive set association and correlation-policy sensitivity, not primary replication or biological validation.',
        source_sha256={p:digest(ROOT/p) for p in sorted(sources)},input_specification_sha256={p:digest(ROOT/p) for p in specs},
        files=[dict(path=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(snap.iterdir()) if p.is_file()])
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Frozen',len(manifest['files']),'files;',len(sources)+len(specs),'source/specification identities')


if __name__=='__main__':main()
