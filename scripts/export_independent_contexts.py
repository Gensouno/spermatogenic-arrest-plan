#!/usr/bin/env python3
"""Freeze this context extension once; earlier snapshots remain immutable."""
import json
import shutil
from check_independent_contexts import TABLES,JSON_FILES,FIGURES,check
from run_independent_contexts import ROOT,digest
from export_reference_limits import local_source_closure


def main():
    out=ROOT/'results/independent_contexts';check(out)
    if json.loads((out/'reproduction.json').read_text())['status']!='passed':raise ValueError('Numerical reproduction required')
    snap=ROOT/'independent_contexts_snapshot'
    if (snap/'manifest.json').exists():raise ValueError('Snapshot already frozen; use a new name')
    snap.mkdir(exist_ok=True)
    for name in TABLES+JSON_FILES+FIGURES+['reproduction.json']:shutil.copy2(out/name,snap/name)
    sources=local_source_closure(ROOT,['scripts/run_independent_contexts.py','scripts/audit_independent_metadata.py',
        'scripts/check_independent_contexts.py','scripts/fetch_independent_contexts.py','scripts/plot_independent_contexts.py',
        'scripts/discover_independent_contexts.py','scripts/export_independent_contexts.py','tests/test_independent_contexts.py'])+['code/analyze_independent_contexts.R']
    specs=['data/independent_discovery_contract.json','data/independent_contexts_contract.json',
           'data/independent_contexts_inputs.json','docs/INDEPENDENT_CONTEXTS.md','envs/python.sh','envs/rscript.sh','envs/network.sh']
    manifest=dict(scope='Four exploratory source contexts and source assertions; no expanded primary family or biological validation.',
        source_sha256={p:digest(ROOT/p) for p in sources},input_specification_sha256={p:digest(ROOT/p) for p in specs},
        files=[dict(path=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(snap.iterdir()) if p.is_file()])
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('Frozen',len(manifest['files']),'files;',len(sources)+len(specs),'source/specification checks')


if __name__=='__main__':main()
