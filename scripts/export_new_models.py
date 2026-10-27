#!/usr/bin/env python3
"""Export bounded new-source diagnostics and source code fingerprints."""
import json
from pathlib import Path
import shutil
from export_reference_limits import local_source_closure
from freeze_replicability_inputs import digest
ROOT=Path(__file__).resolve().parents[1]
TABLES=['sample_evidence.csv','Setx_archive_audit.csv','Eef1g_source_annotation.csv','gene_annotation.csv',
    'cohort_summary.csv','PCA.csv','runtime_R_versions.csv','candidate_context.csv','context_summary.csv','gene_coverage.csv',
    'source_mapping_mass.csv']
PRIVATE_TABLES=['conditional_model_statistics.csv.gz','leave_one_out_effects.csv.gz','TPM_pseudocount_effects.csv.gz']
AUDITS=['preparation_audit.json','design.json']
FIGURES=['fixed_priorities_context.png','sample_QC.png']
SPECS=['data/new_model_audit_contract.json','data/new_model_inputs.json','data/replicability_discovery_contract.json',
       'data/new_model_review.json','envs/rscript.sh','envs/python.sh','envs/network.sh']
ENTRIES=['code/prepare_new_models.py','code/summarize_new_models.py','scripts/run_new_models.py',
    'scripts/check_new_models.py','scripts/check_new_models_reproduction.py','scripts/export_new_models.py',
    'scripts/discover_replicability_data.py','scripts/review_replicability_metadata.py','scripts/fetch_new_model_audit.py',
    'scripts/freeze_new_model_inputs.py','tests/test_new_models.py']


def main():
    snap=ROOT/'new_models_snapshot';snap.mkdir(exist_ok=True);out=ROOT/'results/new_models'
    for n in TABLES+AUDITS+FIGURES:shutil.copy2(out/n,snap/n)
    # Only compact retrieval audit/metadata identifiers, not unbounded SOFT text.
    for n in ['search_inventory.csv','query_audit.json','metadata_review.csv','new_assets.json']:
        shutil.copy2(ROOT/'results/replicability_discovery'/n,snap/n)
    sources=local_source_closure(ROOT,ENTRIES)+['code/analyze_new_models.R']
    manifest={'scope':'Targeted source eligibility audit and exploratory sorted-count / whole-testis TPM contexts; not confirmation of shared mechanism.',
        'source_sha256':{p:digest(ROOT/p) for p in sources},'input_specification_sha256':{p:digest(ROOT/p) for p in SPECS},
        'files':[{'path':str(p.relative_to(snap)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(snap.rglob('*'))
                 if p.is_file() and p.name not in ['manifest.json','.DS_Store']]}
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Exported',len(manifest['files']),'files;',len(sources)+len(SPECS),'source/spec checks')


if __name__=='__main__':main()
