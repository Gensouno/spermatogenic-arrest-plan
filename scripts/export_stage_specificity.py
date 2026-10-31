#!/usr/bin/env python3
"""Freeze this distinct descriptive branch without rewriting preceding evidence."""
import json
import shutil
from check_stage_specificity import TABLES,JSON_FILES,FIGURES,check
from run_stage_specificity import ROOT,digest
from export_reference_limits import local_source_closure


def main():
    out=ROOT/'results/stage_specificity';check(out)
    if json.loads((out/'reproduction.json').read_text())['status']!='passed':raise ValueError('Reproduction required')
    snap=ROOT/'stage_specificity_snapshot'
    if (snap/'manifest.json').exists():raise ValueError('Snapshot frozen; use a new branch name')
    snap.mkdir(exist_ok=True)
    for name in TABLES+JSON_FILES+FIGURES+['reproduction.json']:shutil.copy2(out/name,snap/name)
    sources=local_source_closure(ROOT,['scripts/run_stage_specificity.py','scripts/run_stage_resolution.py','scripts/audit_stage_sources.py',
        'scripts/fetch_stage_sources.py','scripts/check_stage_specificity.py','scripts/plot_stage_specificity.py',
        'scripts/export_stage_specificity.py','tests/test_stage_specificity.py'])
    specs=['data/stage_specificity_contract.json','data/stage_resolution_contract.json','data/stage_specificity_inputs.json',
        'data/stage_source_retrieval.json','docs/STAGE_SPECIFICITY.md','envs/python.sh','envs/network.sh']
    manifest=dict(scope='Post-hoc reference-informed descriptive matching and source metadata reconciliation; not biological validation or an expanded primary family.',
        source_sha256={p:digest(ROOT/p) for p in sources},input_specification_sha256={p:digest(ROOT/p) for p in specs},
        files=[dict(path=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(snap.iterdir()) if p.is_file()])
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('Frozen',len(manifest['files']),'files;',len(sources)+len(specs),'source/specification identities')


if __name__=='__main__':main()
