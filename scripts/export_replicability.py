#!/usr/bin/env python3
"""Compact reproducible replication results; no individual count matrices."""
import json
from pathlib import Path
import shutil
from export_reference_limits import local_source_closure
from freeze_replicability_inputs import digest
ROOT=Path(__file__).resolve().parents[1]
TABLES=['source_statistics_audit.csv','LOO_recomputation_audit.csv','R_adjustment_audit.csv',
    'counts_primary_unit_evidence.csv.gz','legacy_eight_cohort_sensitivity_unit_evidence.csv.gz',
    'directional_replicability_tests.csv.gz','primary_gene_summary.csv.gz','primary_BY_findings.csv',
    'top_ranked_hypotheses.csv','LOO_gene_cohort_summary.csv.gz','LOO_gene_unit_summary.csv.gz',
    'descriptive_1p5fold_candidates.csv.gz','direction_patterns.csv','replicability_summary.csv',
    'descriptive_support_counts.csv','gene_coverage.csv','all_primary_gene_stage_context.csv.gz',
    'candidate_stage_context.csv','stable_descriptive_priorities.csv','stage_context_summary.csv',
    'candidate_unit_support.csv','stage_source_cross_counts.csv','runtime_R_versions.csv',
    'external_candidate_context.csv','external_context_summary.csv','external_context_coverage.csv']
AUDITS=['design.json','stage_annotation_design.json','runtime_Python_versions.json','external_context_design.json']
FIGURES=['descriptive_support.png','stable_priority_context.png','external_context_coverage.png']
PRIVATE_TABLES=['recomputed_study_statistics.csv.gz','recomputed_count_LOO_effects.csv.gz']
ENTRIES=['code/analyze_replicability.py','code/annotate_replicability.py','code/plot_replicability.py','code/check_replicability_context.py',
    'scripts/freeze_replicability_inputs.py','scripts/run_replicability.py','scripts/export_replicability.py',
    'scripts/check_replicability.py','scripts/check_replicability_reproduction.py','tests/test_replicability.py']
R_SOURCES=['code/audit_replicability_statistics.R','code/audit_replicability_LOO.R','code/audit_replicability_adjustments.R']
SPECS=['data/replicability_contract.json','data/replicability_inputs.json','data/replicability_context_contract.json',
       'envs/rscript.sh','envs/python.sh','envs/network.sh']


def main():
    out=ROOT/'results/replicability';snap=ROOT/'replicability_snapshot';snap.mkdir(exist_ok=True)
    for n in TABLES+AUDITS+FIGURES:shutil.copy2(out/n,snap/n)
    sources=local_source_closure(ROOT,ENTRIES)+R_SOURCES
    manifest={'scope':'Directional hypothesis replication and descriptive stage-context priorities, not a validated mechanism.',
        'source_sha256':{p:digest(ROOT/p) for p in sources},'input_specification_sha256':{p:digest(ROOT/p) for p in SPECS},
        'files':[{'path':str(p.relative_to(snap)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(snap.rglob('*'))
                 if p.is_file() and p.name not in ['manifest.json','.DS_Store']]}
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Exported',len(manifest['files']),'files;',len(sources)+len(SPECS),'source/spec checks')


if __name__=='__main__':main()
