#!/usr/bin/env python3
"""Export compact results and bind the exact source/specification files."""
import json
from pathlib import Path
import shutil

from emptydrops_range_client import digest

ROOT=Path(__file__).resolve().parents[1]
TABLES=['profile_evidence.csv','group_inventory.csv','cell_complexity_summary.csv','holdout_metrics.csv',
    'full_fit_weights.csv','fold_weights.csv','stratum_summary.csv','early_meiotic_comparison.csv',
    'matched_count_branch_comparison.csv','design_conditions.csv']
AUDITS=['preparation_audit.json','design.json']
SOURCES=['scripts/emptydrops_range_client.py','scripts/fetch_emptydrops_p15.py','scripts/run_emptydrops_p15.py',
    'scripts/check_emptydrops_p15.py','scripts/check_emptydrops_reproduction.py','scripts/export_emptydrops_p15.py',
    'code/prepare_emptydrops_p15.py','code/diagnose_emptydrops_p15.py','tests/test_emptydrops_p15.py',
    'code/prepare_labeled_reference.py','code/diagnose_reference_transfer.py','code/check_sorted_expansion.py',
    'scripts/discover_expansion.py','scripts/fetch_labeled_reference.py','scripts/check_project.py',
    'data/emptydrops_p15_contract.json']
SPECS=['data/emptydrops_p15_assets.json','data/labeled_reference_assets.json','labeled_reference_snapshot/design.json']


def main():
    out=ROOT/'results/emptydrops_p15';snap=ROOT/'emptydrops_p15_snapshot';snap.mkdir(exist_ok=True)
    for name in TABLES+AUDITS+['gene_design.csv.gz','early_meiotic_diagnostics.png']:shutil.copy2(out/name,snap/name)
    manifest={'scope':json.loads((out/'design.json').read_text())['scope'],
        'source_sha256':{p:digest(ROOT/p) for p in SOURCES},
        'input_specification_sha256':{p:digest(ROOT/p) for p in SPECS},
        'files':[{'path':str(p.relative_to(snap)),'bytes':p.stat().st_size,'sha256':digest(p)}
            for p in sorted(snap.rglob('*')) if p.is_file() and p.name not in ['manifest.json','.DS_Store']]}
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Exported complete-P15 snapshot:',len(manifest['files']),'files')


if __name__=='__main__':main()
