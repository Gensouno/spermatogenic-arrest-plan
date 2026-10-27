#!/usr/bin/env python3
"""Compact independent-source diagnostic snapshot and exact source fingerprints."""
import hashlib
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[1]
TABLES=['sample_evidence.csv','candidate_review.csv','replicate_concordance.csv','fold_weights.csv','full_fit_weights.csv',
    'holdout_metrics.csv','stage_comparison.csv','design_conditions.csv','gene_design.csv.gz']
AUDITS=['mapping_audit.json','design.json']
SOURCES=['code/diagnose_leptotene.py','scripts/discover_leptotene.py','scripts/fetch_leptotene.py','scripts/run_leptotene.py',
    'scripts/export_leptotene.py','scripts/check_leptotene.py','scripts/check_leptotene_reproduction.py','tests/test_leptotene.py',
    'code/diagnose_reference_transfer.py','code/check_sorted_expansion.py','code/prepare_expansion.py',
    'scripts/discover_expansion.py','scripts/check_project.py','data/leptotene_contract.json']
SPECS=['data/leptotene_assets.json','data/leptotene_paper_asset.json','data/labeled_reference_assets.json','labeled_reference_snapshot/design.json']


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    out=ROOT/'results/leptotene_followup';snap=ROOT/'leptotene_snapshot';snap.mkdir(exist_ok=True)
    for name in TABLES+AUDITS+['stage_diagnostics.png']:shutil.copy2(out/name,snap/name)
    manifest={'scope':json.loads((out/'design.json').read_text())['scope'],
        'source_sha256':{p:digest(ROOT/p) for p in SOURCES},'input_specification_sha256':{p:digest(ROOT/p) for p in SPECS},
        'files':[{'path':str(p.relative_to(snap)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(snap.rglob('*')) if p.is_file() and p.name not in ['manifest.json','.DS_Store']]}
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Exported independent-source snapshot:',len(manifest['files']),'files',flush=True)


if __name__=='__main__':main()
