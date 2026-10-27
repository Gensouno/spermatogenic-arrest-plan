#!/usr/bin/env python3
"""Publish compact diagnostics and recursively fingerprint local Python sources."""
import ast
import json
from pathlib import Path
import shutil

from freeze_reference_limits_inputs import digest

ROOT=Path(__file__).resolve().parents[1]
TABLES=['chen_profile_evidence.csv','chen_cell_QC.csv.gz','chen_library_inventory.csv','chen_stage_inventory.csv','chen_batch_evidence.csv','chen_unresolved_cells.csv',
    'gene_union_exon_lengths.csv.gz','technical_full_fit_weights.csv','technical_holdout_metrics.csv','technical_comparison.csv',
    'technical_conditions.csv','technical_gene_design.csv.gz','technical_fold_weights.csv.gz','target_evidence.csv',
    'training_marker_panel.csv.gz','transfer_full_weights.csv','transfer_fold_weights.csv.gz','transfer_holdout_metrics.csv',
    'transfer_comparison.csv','transfer_conditions.csv','prototype_evidence.csv','transfer_label_summary.csv','transfer_gene_folds.csv.gz',
    'source_library_holdout.csv','source_library_holdout_folds.csv.gz','source_library_conditions.csv',
    'source_library_training_inventory.csv','source_library_external_ranges.csv','technical_focal_summary.csv',
    'transfer_focal_summary.csv','source_holdout_controls.csv','legacy_sorted_comparison.csv','legacy_sorted_fold_weights.csv.gz',
    'legacy_sorted_conditions.csv','legacy_sorted_samples.csv','legacy_sorted_gene_design.csv.gz']
AUDITS=['chen_preparation_audit.json','technical_design.json','transfer_design.json','source_library_design.json','decision_summary.json','legacy_sorted_design.json']
ENTRY_SOURCES=['scripts/fetch_chen_umi.py','scripts/run_reference_limits.py','scripts/freeze_reference_limits_inputs.py',
    'scripts/export_reference_limits.py','scripts/check_reference_limits.py','scripts/check_reference_limits_reproduction.py',
    'code/prepare_chen_reference.py','code/diagnose_reference_limits.py','code/check_chen_transfer.py','code/check_chen_library_holdout.py',
    'code/plot_reference_limits.py','code/check_legacy_sorted_transfer.py','tests/test_reference_limits.py']
SPECS=['data/reference_limits_contract.json','data/reference_limits_inputs.json','data/reference_limits_sorted_inputs.json','data/chen_umi_assets.json','data/leptotene_assets.json']


def local_source_closure(root,entries):
    found=set();pending=list(entries)
    while pending:
        relative=pending.pop()
        if relative in found:continue
        path=root/relative
        if not path.is_file():raise FileNotFoundError('Missing source '+relative)
        found.add(relative)
        tree=ast.parse(path.read_text())
        modules=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.ImportFrom) and node.module:modules.append(node.module.split('.')[0])
            if isinstance(node,ast.Import):modules.extend(a.name.split('.')[0] for a in node.names)
        for module in modules:
            matches=[str(p.relative_to(root)) for folder in ['code','scripts','tests'] if (p:=root/folder/(module+'.py')).is_file()]
            if len(matches)>1:raise ValueError('Ambiguous local module '+module)
            pending+=matches
    return sorted(found)


def main():
    out=ROOT/'results/reference_limits';snap=ROOT/'reference_limits_snapshot';snap.mkdir(exist_ok=True)
    for name in TABLES+AUDITS+['reference_limit_diagnostics.png','candidate_transfer_diagnostics.png']:shutil.copy2(out/name,snap/name)
    sources=local_source_closure(ROOT,ENTRY_SOURCES)
    manifest={'scope':'Technical sensitivities, reconciled Chen UMI resource and exploratory transfer failure controls; not a calibrated cell-fraction method.',
        'source_sha256':{p:digest(ROOT/p) for p in sources},'input_specification_sha256':{p:digest(ROOT/p) for p in SPECS},
        'files':[{'path':str(p.relative_to(snap)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(snap.rglob('*'))
            if p.is_file() and p.name not in ['manifest.json','.DS_Store']]}
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Exported',len(manifest['files']),'files;',len(sources),'local Python source fingerprints',flush=True)


if __name__=='__main__':main()
