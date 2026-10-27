#!/usr/bin/env python3
"""Freeze the new branch only; never alter the older snapshot families."""
import json
from pathlib import Path
import shutil
from check_setx_interaction import TABLES,AUDITS,FIGURES,check_results
from export_reference_limits import local_source_closure
from run_setx_interaction import digest

ROOT=Path(__file__).resolve().parents[1]


def main():
    check_results(ROOT)
    snap=ROOT/'setx_interaction_snapshot'
    if (snap/'manifest.json').exists():
        raise ValueError('Already frozen; create an explicitly new snapshot for subsequent changes')
    snap.mkdir(exist_ok=True)
    out=ROOT/'results/setx_interaction'
    for name in TABLES+AUDITS+FIGURES+['reproduction.json']:
        shutil.copy2(out/name,snap/name)
    sources=local_source_closure(ROOT,['scripts/run_setx_interaction.py','scripts/check_setx_interaction.py',
        'scripts/export_setx_interaction.py','tests/test_setx_interaction.py'])+['code/analyze_setx_interaction.R']
    specs=['data/setx_interaction_contract.json','data/setx_interaction_inputs.json',
        'envs/python.sh','envs/rscript.sh','envs/network.sh','docs/SETX_INTERACTION.md']
    manifest=dict(scope='Post-hoc conditional source-library age-by-genotype interaction; not confirmatory biological validation.',
        source_sha256={p:digest(ROOT/p) for p in sources},
        input_specification_sha256={p:digest(ROOT/p) for p in specs},
        files=[dict(path=str(p.relative_to(snap)),bytes=p.stat().st_size,sha256=digest(p))
            for p in sorted(snap.iterdir()) if p.is_file() and p.name!='.DS_Store'])
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Frozen',len(manifest['files']),'files;',len(sources)+len(specs),'source/specification hashes')


if __name__=='__main__':main()
