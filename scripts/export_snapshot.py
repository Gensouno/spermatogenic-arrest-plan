#!/usr/bin/env python3
"""Export a small, explicit, inspectable result set; no raw matrices or logs."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results_snapshot'
ASSETS = [
    'all_cohort_comparison/common_gene_spearman.csv',
    'all_cohort_comparison/pairwise_gene_spearman.csv',
    'all_cohort_comparison/pairwise_gene_numbers.csv',
    'all_cohort_comparison/diagnostic_markers.csv',
    'all_cohort_comparison/eight_cohort_effects.png',
    'all_cohort_comparison/within_study_PCA.png',
    'all_cohort_comparison/PCA_sample_labels.csv',
    'import_qc/official_sample_evidence.csv', 'import_qc/metadata_sources.json',
    'import_qc/count_cohorts.csv', 'import_qc/extended_rnaseq_cohorts.csv',
    'import_qc/microarray_cohorts.csv', 'import_qc/GSE101511_mapping_mass.csv',
    'import_qc/R_sessionInfo.txt', 'import_qc/extended_R_sessionInfo.txt',
    'import_qc/microarray_R_sessionInfo.txt', 'robustness/R_sessionInfo.txt',
    'robustness/baseline_reproduction.csv', 'robustness/sensitivity_metrics.csv',
    'robustness/sensitivity_summary.csv', 'robustness/pair_ranges.csv',
    'robustness/pairwise_sensitivity.csv', 'robustness/key_pair_sensitivity.png',
    'reference_calibration/summary.csv', 'reference_calibration/per_batch.csv',
    'reference_calibration/design.json', 'reference_calibration/official_atlas_sources.csv',
    'reference_calibration/atlas_metadata_audit.json',
    'reference_calibration/cross_source_error.png',
    'mixture_audit/calibration_summary.csv', 'mixture_audit/real_heldout.csv',
    'mixture_audit/provenance.json', 'mixture_audit/mixture_audit.png',
    'composition_prediction/prediction_metrics.csv',
    'composition_prediction/covariance_terms.csv',
    'external_sorted_reference/fitted_RNA_weights.csv',
    'external_sorted_reference/within_meiotic_stage_contrast.csv',
    'external_sorted_reference/excluded_nonnumeric_rows.csv',
    'external_sorted_reference/provenance.json',
    'stage_reference/fitted_weights.csv', 'stage_reference/crossfit_metrics.csv',
    'stage_reference/excluded_excel_date_gene_rows.csv',
    'stage_reference/germ_cluster_sizes.csv', 'stage_reference/provenance.json',
    'stage_reference/controlled_resolution_comparison.png',
]
LABELS = {'GSE236388': 'XAP5', 'GSE283904': 'N6amt1', 'GSE68283': 'Rfx2_P29_30',
          'GSE68284': 'Rfx2_P21', 'GSE101511': 'Hnrnpgt', 'GSE294209': 'Rnu11',
          'GSE38096': 'Sin3a', 'GSE42601': 'Dazap1'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_snapshot():
    manifest = json.loads((OUT / 'manifest.json').read_text())
    declared = {row['path'] for row in manifest['files']}
    actual = {str(p.relative_to(OUT)) for p in OUT.rglob('*') if p.is_file() and p.name != '.DS_Store'} - {'manifest.json'}
    if declared != actual:
        raise ValueError(f'Snapshot inventory differs: {declared ^ actual}')
    for row in manifest['files']:
        path = (OUT / row['path']).resolve()
        if not path.is_relative_to(OUT.resolve()):
            raise ValueError('Snapshot path escapes output directory')
        if path.stat().st_size != row['bytes'] or digest(path) != row['sha256']:
            raise ValueError(f'Snapshot checksum failed: {row["path"]}')
    for relative, checksum in manifest['source_sha256'].items():
        path = (ROOT / relative).resolve()
        if not path.is_relative_to(ROOT.resolve()) or digest(path) != checksum:
            raise ValueError(f'Snapshot source fingerprint differs: {relative}; review and regenerate')
    for relative, key in [('data/analysis_assets.json', 'asset_manifest_sha256'),
                          ('data/samples.tsv', 'sample_sheet_sha256')]:
        if digest(ROOT / relative) != manifest[key]:
            raise ValueError(f'Snapshot input specification differs: {relative}')
    print('Snapshot integrity OK:', len(declared), 'files,', sum(r['bytes'] for r in manifest['files']), 'bytes')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.verify:
        verify_snapshot()
        return
    clean = json.loads((ROOT / 'results/reproducibility/clean_run.json').read_text())
    if clean['status'] != 'complete' or not all(r['matches'] for r in clean['compared_tables']):
        raise ValueError('Successful clean-directory reproduction is required before export')
    completed = []
    for path in (ROOT / 'logs/analysis').glob('*/run.json'):
        record = json.loads(path.read_text())
        if (record.get('status') == 'complete' and len(record['stages']) == 19
                and not record['asset_check_skipped']):
            completed.append((path, record))
    if not completed:
        raise ValueError('A complete, asset-checked 19-stage baseline run is required')
    run_path, run = sorted(completed, key=lambda item: str(item[0]))[-1]
    # Packaging helpers may evolve without altering numerical analyses. Bind
    # results to every analysis source and the actual runner/input checker.
    analysis_sources = {p: h for p, h in run['source_sha256'].items()
                        if p.startswith('code/') or p in ['scripts/run_analysis.py', 'scripts/data_assets.py']}
    for relative, checksum in {**analysis_sources, **run['curation_sha256']}.items():
        if digest(ROOT / relative) != checksum:
            raise ValueError(f'Analysis changed since the completed run: {relative}; rerun before export')
    OUT.mkdir(exist_ok=True)
    for relative in ASSETS:
        source = ROOT / 'results' / relative
        target = OUT / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if source.suffix == '.txt':
            # Session reports use padded columns; preserve their information
            # while keeping generated snapshot text free of trailing spaces.
            target.write_text('\n'.join(line.rstrip() for line in source.read_text().splitlines()) + '\n')
        else:
            shutil.copyfile(source, target)
    counts = pd.read_csv(ROOT / 'results/import_qc/count_cohorts.csv').assign(unit='counts')
    extended = pd.read_csv(ROOT / 'results/import_qc/extended_rnaseq_cohorts.csv')
    arrays = pd.read_csv(ROOT / 'results/import_qc/microarray_cohorts.csv').assign(unit='microarray')
    columns = ['gse', 'unit', 'WT', 'KO', 'tested_genes', 'DEG_FDR05_absLFC1']
    summary = pd.concat([counts[columns], extended[columns], arrays[columns]], ignore_index=True)
    summary.insert(1, 'perturbation', summary.gse.map(LABELS))
    summary.to_csv(OUT / 'cohort_summary.csv', index=False)
    evidence = pd.read_csv(ROOT / 'results/import_qc/official_sample_evidence.csv')
    if len(evidence) != 56 or evidence.gsm.nunique() != 56 or evidence.gse.nunique() != 8:
        raise ValueError('Unexpected primary sample inventory')
    report = {'scope': clean['scope'], 'status': clean['status'],
              'elapsed_seconds': clean['elapsed_seconds'], 'input_files_copied': clean['input_files_copied'],
              'compared_tables': clean['compared_tables'],
              'primary_samples': len(evidence), 'datasets': evidence.gse.nunique(),
              'genetic_perturbations': 7,
              'baseline_run_id': run_path.parent.name,
              'analysis_source_sha256': analysis_sources,
              'python_version': run['python_version'],
              'cell_number_fractions_validated': False,
              'causal_composition_percentages_validated': False,
              'warning': 'Successful reproducibility checks do not establish biological validity. Same installed dependencies were reused.'}
    (OUT / 'validation_summary.json').write_text(json.dumps(report, indent=2) + '\n')
    files = []
    expected = {'README.md', 'cohort_summary.csv', 'validation_summary.json', *ASSETS}
    present = {str(p.relative_to(OUT)) for p in OUT.rglob('*') if p.is_file() and p.name not in ['manifest.json', '.DS_Store']}
    if present != expected:
        raise ValueError(f'Review unexpected or missing snapshot files: {present ^ expected}')
    for relative in sorted(expected):
        path = OUT / relative
        files.append({'path': relative, 'bytes': path.stat().st_size, 'sha256': digest(path)})
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    sources = {}
    for folder in ['code', 'scripts', 'tests', 'envs']:
        for path in (ROOT / folder).glob('*'):
            if path.is_file() and path.suffix in ['.py', '.R', '.sh', '.txt']:
                sources[str(path.relative_to(ROOT))] = digest(path)
    manifest = {'schema': 1, 'exported_utc': datetime.now(timezone.utc).isoformat(),
                'git_base_commit': commit,
                'git_note': 'Base commit is provenance only; source_sha256 identifies the prepared working-tree code.',
                'asset_manifest_sha256': digest(ROOT / 'data/analysis_assets.json'),
                'sample_sheet_sha256': digest(ROOT / 'data/samples.tsv'),
                'source_sha256': sources, 'files': files}
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    verify_snapshot()


if __name__ == '__main__':
    main()
