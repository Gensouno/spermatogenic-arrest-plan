#!/usr/bin/env python3
"""Compare optional expansion with a frozen baseline; export a separate snapshot."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/data_expansion'
DECISIONS = {
 'GSE259242': ('analyzed_new_arrest_comparison', 'P9 Dcaf2; explicit count columns only; 3 control/3 cKO.'),
 'GSE117706': ('analyzed_HET_control_sensitivity', 'P2 and P7 separately, each 3 control/3 KO; controls are Pramef12+/-, never renamed WT.'),
 'GSE235777': ('analyzed_phenotype_negative_comparator', 'P45 Msl3 3 WT/3 cKO, STAR unstranded column 2; not a guaranteed transcriptional-null truth set.'),
 'GSE232415': ('whole_testis_blocked_sorted_WT_only', '15 internal-label mismatches and 6 whole-testis files in duplicate groups; no inferred relabeling. Four separately checked sorted WT PS/RS profiles used only for reference diagnostics.'),
 'GSE158471': ('analyzed_sorted_WT_reference_check', 'Four Cnot4 WT LZ/PD profiles, FPKM not counts; mm9 symbols and two replicates per stage.'),
 'GSE290629': ('analyzed_sorted_WT_reference_check', 'Twelve author-labeled WT stage profiles; no assertion that stages are independent animals.'),
 'GSE278121': ('unit_pending', 'About 48% matrix entries are non-integer despite raw_counts filename; quantification scale unresolved; not rounded.'),
 'GSE290783': ('unit_and_sample_metadata_pending', 'About 95% entries non-integer; genotype characteristic is literally tissue for all samples; not admitted to primary counts analysis.'),
 'GSE121904': ('age_reference_candidate', 'PND6/14/18/25/30 and adult; capture preparation differs and adult enrichment must be respected; cell-label mapping pending.'),
 'GSE230746': ('age_reference_candidate_resource_deferred', 'P15 WT has two biological replicates and author-labeled Seurat object, but object is 1,122,924,799 bytes compressed; not downloaded.'),
 'GSE228722': ('replicated_reference_candidate', 'P24 scRNA-seq, 4 control/2 NELF-B cKO/3 TDP-43 cKO; not bulk biological replication by summing cells.'),
 'GSE180204': ('priority_future_bulk_candidate', 'P1 Nanos3-Cre and P3 Dhh-Cre arms, 4 vs 4 each; different targeted cell types and floxed controls, not interchangeable.'),
 'GSE190646': ('FPKM_sensitivity_candidate', 'P3/P6/P8 whole-testis FPKM with small replicate counts plus sorted/single-cell components; not primary integer counts.'),
 'GSE172171': ('synchronized_stage_candidate', 'WT/KO early and late pachytene after WIN18,446/RA synchronization; not untreated whole testis.'),
 'GSE230743': ('exclude_unmanipulated_whole_testis_pool', 'Developmentally synchronized preleptotene-enriched testes.'),
 'GSE230747': ('superseries_not_independent_dataset', 'Parent of GSE230743 and GSE230746; do not double count.'),
 'GSE206018': ('exclude_whole_testis_pool', 'Primary Sertoli cells, not whole testis.'),
 'GSE217685': ('exclude_whole_testis_pool', 'THY1-positive sorted spermatogonia.'),
 'GSE252121': ('exclude_whole_testis_pool', 'Sorted Thy1+/cKit+ RNA-seq and ChIP-seq mixed in the series.'),
 'GSE264025': ('exclude_whole_testis_pool', 'Prepared germ cells from P15, not unselected whole testis.'),
 'GSE95225': ('exclude_non_expression_assay', 'Official sample metadata identifies H3K27me3 ChIP/input, despite expression-like study title.'),
 'GSE341157': ('exclude_non_testis_expression_assay', 'Official series is ChIL-seq in NIH3T3 cells; misleading search snippets do not establish whole-testis RNA-seq.'),
 'GSE232398': ('exclude_non_expression_assay', 'ChIP-seq binding profiles.'),
 'GSE232416': ('exclude_non_expression_assay', 'eCLIP-seq RNA binding, not bulk gene expression.'),
 'GSE232417': ('exclude_non_expression_assay', 'ssDRIP-seq R-loop profiles.'),
 'GSE182518': ('mixed_assay_sorted_candidate', 'Sorted CB/RS/spermatocyte and small-RNA/TE/gene outputs; not added as whole testis.'),
 'GSE181132': ('small_sample_candidate', 'Adult whole testis 2 WT/2 Ppp2ca cKO; matrix units still require file inspection.'),
 'GSE242235': ('mixed_assay_candidate', 'Whole-testis RNA arms exist alongside Capture-HiC/ChIP and homozygous/heterozygous genotypes; expression file availability not yet fully checked.'),
 'GSE63166': ('normalized_sensitivity_candidate', 'Multiple small-RNA/expression experiments and processed CPM tables; not verified raw counts.'),
 'GSE79106': ('expression_asset_pending', 'Whole-series supplement list is empty; sample-level files and raw availability still need inspection.'),
 'GSE117707': ('single_cell_candidate', 'PRAMEF12 scRNA-seq, not extra bulk replicates; source labels still need extraction.'),
 'GSE151564': ('reference_candidate_pending', 'AR-signaling study; not admitted without checking population, genotype and expression format.'),
}


def main():
    baseline = pd.read_csv(ROOT / 'results/all_cohort_comparison/common_logFC.csv.gz', index_col=0)
    summary = pd.read_csv(OUT / 'cohort_summary.csv')
    effects = {c: pd.read_csv(OUT / f'de/{c}_de.csv').set_index('gene').logFC for c in summary.cohort}
    pairs = []
    for c, effect in effects.items():
        common = baseline.index.intersection(effect.index)
        for old in baseline:
            pairs.append({'expansion_cohort': c, 'baseline_cohort': old, 'genes': len(common),
                          'spearman': effect.loc[common].corr(baseline.loc[common, old], method='spearman')})
    pd.DataFrame(pairs).to_csv(OUT / 'fixed_core_gene_correlations.csv', index=False)
    loo = pd.read_csv(OUT / 'leave_one_out_logFC.csv.gz')
    sensitivity = []
    for (cohort, deleted), frame in loo.groupby(['cohort', 'deleted']):
        effect = frame.set_index('gene').logFC
        common = baseline.index.intersection(effect.index).intersection(effects[cohort].index)
        for old in baseline:
            sensitivity.append({'expansion_cohort': cohort, 'deleted': deleted, 'baseline_cohort': old,
                'genes': len(common), 'baseline_rho_on_same_genes': effects[cohort].loc[common].corr(baseline.loc[common, old], method='spearman'),
                'leave_one_out_rho': effect.loc[common].corr(baseline.loc[common, old], method='spearman')})
    sensitivity = pd.DataFrame(sensitivity)
    sensitivity.to_csv(OUT / 'cross_cohort_leave_one_out.csv', index=False)
    sensitivity.groupby(['expansion_cohort', 'baseline_cohort']).leave_one_out_rho.agg(['min', 'max']).to_csv(OUT / 'cross_cohort_ranges.csv')
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    frame = pd.DataFrame(pairs).pivot(index='expansion_cohort', columns='baseline_cohort', values='spearman').reindex(summary.cohort)
    labels = ['Dcaf2 P9', 'Pramef12 P2 [HET ctrl]', 'Pramef12 P7 [HET ctrl]', 'Msl3 P45 [pheno-negative]']
    fig, ax = plt.subplots(figsize=(10, 4.8), layout='constrained')
    im = ax.imshow(frame, cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set(xticks=range(len(frame.columns)), xticklabels=frame.columns,
           yticks=range(4), yticklabels=labels, title='Optional expansion vs frozen baseline gene universe')
    ax.tick_params(axis='x', rotation=55)
    for i in range(len(frame)):
        for j in range(len(frame.columns)):
            ax.text(j, i, f'{frame.iloc[i,j]:.2f}', ha='center', va='center', fontsize=9)
    fig.colorbar(im, ax=ax, shrink=.8, label='Descriptive Spearman')
    fig.savefig(OUT / 'expansion_effects.png', dpi=180)
    plt.close(fig)
    records = []
    for path in sorted(OUT.glob('GSE*_metadata.json')):
        record = json.loads(path.read_text()); gse = record['accession']
        decision, reason = DECISIONS.get(gse, ('unreviewed', 'Not enough evidence for inclusion'))
        records.append({'gse': gse, 'title': record['series']['title'], 'samples_in_series': len(record['samples']),
            'decision': decision, 'reason': reason, 'metadata_url': record['source_url'], 'metadata_sha256': record['sha256']})
    review = {'review_date_local': '2026-10-23', 'scope': 'Targeted search, not an exhaustive systematic review or preregistration.',
              'selection_basis': 'Source tissue, controls, assay, units and file identity; no model selected on correlation with baseline.', 'candidates': records}
    (ROOT / 'data/expansion_review.json').write_text(json.dumps(review, indent=2) + '\n')
    snap = ROOT / 'expansion_snapshot'; snap.mkdir(exist_ok=True)
    selected = ['cohort_summary.csv', 'sample_sensitivity.csv', 'fixed_core_gene_correlations.csv',
        'cross_cohort_ranges.csv', 'GSE232415_file_identity_audit.csv', 'unresolved_count_units.csv',
        'expansion_effects.png', 'sorted_validation/summary.csv', 'sorted_validation/fitted_weights.csv',
        'sorted_validation/source_sample_evidence.csv', 'sorted_validation/design.json',
        'sorted_validation/new_sorted_compatibility.png']
    for relative in selected:
        target = snap / relative; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(OUT / relative, target)
    for c in summary.cohort:
        shutil.copyfile(OUT / f'{c}_sample_evidence.csv', snap / f'{c}_sample_evidence.csv')
    sources = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in [ROOT / 'code' / name for name in ['prepare_expansion.py', 'de_expansion.R',
                   'map_sorted_expansion.R', 'check_sorted_expansion.py', 'summarize_expansion.py']]}
    for name in ['scripts/discover_expansion.py', 'scripts/fetch_expansion_asset.py',
                 'scripts/run_expansion.py', 'scripts/check_expansion_reproduction.py',
                 'tests/test_expansion.py', 'tests/test_proxy_env.py',
                 'envs/network.sh', 'envs/python.sh', 'envs/rscript.sh']:
        sources[name] = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
    inventory = [{'path': str(p.relative_to(snap)), 'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                 for p in sorted(snap.rglob('*')) if p.is_file() and p.name != 'manifest.json']
    (snap / 'manifest.json').write_text(json.dumps({'generated_utc': datetime.now(timezone.utc).isoformat(),
        'baseline': 'Frozen eight-cohort analysis, commit 7d8e45e; not overwritten by expansion.',
        'baseline_effect_matrix_sha256': hashlib.sha256((ROOT / 'results/all_cohort_comparison/common_logFC.csv.gz').read_bytes()).hexdigest(),
        'whole_testis_samples_added': 24, 'bulk_contrasts_added': 4, 'sorted_WT_profiles_checked': 20,
        'source_sha256': sources, 'asset_manifest_sha256': hashlib.sha256((ROOT / 'data/expansion_assets.json').read_bytes()).hexdigest(),
        'review_sha256': hashlib.sha256((ROOT / 'data/expansion_review.json').read_bytes()).hexdigest(),
        'files': inventory, 'cell_proportions_validated': False}, indent=2) + '\n')
    print(frame.to_string())
    print('Candidates reviewed:', len(records), '; additional whole-testis samples: 24; sorted WT profiles: 20')


if __name__ == '__main__':
    main()
