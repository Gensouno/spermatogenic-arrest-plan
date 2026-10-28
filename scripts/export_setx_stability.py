#!/usr/bin/env python3
"""Export one new compact snapshot without altering preceding branches."""
import hashlib
import json
from pathlib import Path
import shutil
from check_setx_stability import TABLES,JSON_FILES,FIGURES,check
from export_reference_limits import local_source_closure

ROOT=Path(__file__).resolve().parents[1]


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    check(ROOT);snap=ROOT/'setx_stability_snapshot'
    if (snap/'manifest.json').exists():raise ValueError('Snapshot already frozen; use a new name for subsequent changes')
    snap.mkdir(exist_ok=True)
    out=ROOT/'results/setx_stability'
    for name in TABLES+JSON_FILES+FIGURES+['reproduction.json']:shutil.copy2(out/name,snap/name)
    sources=local_source_closure(ROOT,['scripts/run_setx_stability.py','scripts/check_setx_stability.py',
        'scripts/export_setx_stability.py','tests/test_setx_stability.py'])+['code/analyze_setx_stability.R','code/analyze_setx_length.R','code/analyze_setx_interaction.R']
    specs=['data/setx_stability_contract.json','data/setx_stability_inputs.json',
        'envs/python.sh','envs/rscript.sh','envs/network.sh','docs/SETX_STABILITY.md']
    manifest=dict(scope='Same-source post-hoc omission stability and matched-universe diagnostics, not biological validation.',
        source_sha256={p:digest(ROOT/p) for p in sources},input_specification_sha256={p:digest(ROOT/p) for p in specs},
        files=[dict(path=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(snap.iterdir()) if p.is_file()])
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Frozen',len(manifest['files']),'files;',len(sources)+len(specs),'source/specification checks')


if __name__=='__main__':main()
