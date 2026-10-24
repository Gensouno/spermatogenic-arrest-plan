#!/usr/bin/env python3
"""Export bounded Lrh1 diagnostics without changing either prior snapshot."""
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/followup'
COHORTS = ['Lrh1_P1_NanosCre', 'Lrh1_P3_DhhCre']


def unique_symbol_effects(frame):
    mapped = frame.dropna(subset=['symbol']).copy()
    # Exclude all ambiguous symbols, not the gene with the less convenient effect.
    mapped = mapped[~mapped.symbol.duplicated(keep=False)]
    return mapped.set_index('symbol').logFC


def main():
    baseline = pd.read_csv(ROOT / 'results/all_cohort_comparison/common_logFC.csv.gz', index_col=0)
    effects = {c: pd.read_csv(OUT / f'de/{c}.csv') for c in COHORTS}
    mapping = pd.read_csv(OUT / 'Lrh1_mapping.csv')
    gene_mapping = dict(zip(mapping.input_id.str.replace(r'\.[0-9]+$', '', regex=True), mapping.symbol))
    correlations, sensitivity = [], []
    loo = pd.read_csv(OUT / 'leave_one_out_effects.csv.gz')
    for cohort, de in effects.items():
        effect = unique_symbol_effects(de)
        genes = baseline.index.intersection(effect.index)
        for old in baseline:
            correlations.append({'cohort': cohort, 'baseline_cohort': old, 'genes': len(genes),
                'spearman': effect.loc[genes].corr(baseline.loc[genes, old], method='spearman')})
        for deleted, part in loo[loo.cohort == cohort].groupby('deleted'):
            part = part.assign(symbol=part.gene.map(gene_mapping))
            refit = unique_symbol_effects(part)
            shared = genes.intersection(refit.index)
            for old in baseline:
                sensitivity.append({'cohort': cohort, 'deleted': deleted, 'baseline_cohort': old,
                    'genes': len(shared),
                    'full_rho_same_genes': effect.loc[shared].corr(baseline.loc[shared, old], method='spearman'),
                    'leave_one_out_rho': refit.loc[shared].corr(baseline.loc[shared, old], method='spearman')})
    corr = pd.DataFrame(correlations);corr.to_csv(OUT / 'baseline_correlations.csv', index=False)
    sensitivity = pd.DataFrame(sensitivity)
    sensitivity.to_csv(OUT / 'baseline_leave_one_out.csv', index=False)
    sensitivity.groupby(['cohort', 'baseline_cohort']).leave_one_out_rho.agg(['min', 'max']).to_csv(OUT / 'baseline_ranges.csv')
    workbook = ROOT / 'data/raw/followup/GSE180204/GSE180204_Agrimson_Lrh1_Supplemental_Tables.xlsx'
    author = pd.read_excel(workbook, sheet_name='Nanos P1 DEGs')
    contract = json.loads((ROOT / 'data/followup_contract.json').read_text())
    targets = contract['targets']['GSE180204']['paper_direction_checks']
    if set(author['Gene Symbol']) != set(targets):
        raise ValueError('Paper P1 gene table differs from the prespecified four-gene check')
    primary = effects[COHORTS[0]].set_index('gene')
    alternate = pd.read_csv(OUT / f'de/{COHORTS[0]}_B2_only.csv').set_index('gene')
    qcheck = pd.read_csv(OUT / f'de/{COHORTS[0]}_edgeRQL.csv').set_index('gene')
    checks = []
    for _, r in author.iterrows():
        gene = r['Gene ID']; fit = primary.loc[gene]
        expected = 1 if targets[r['Gene Symbol']] == 'up' else -1
        checks.append({'gene_id': gene, 'author_symbol': r['Gene Symbol'],
            'author_logFC': r['Log2 Fold Change'], 'author_adjusted_P': r['BH adjusted p-value'],
            'reanalyzed_logFC': fit.logFC, 'reanalyzed_adjusted_P': fit['adj.P.Val'],
            'B2_only_logFC': alternate.loc[gene, 'logFC'],
            'edgeRQL_logFC': qcheck.loc[gene, 'logFC'], 'edgeRQL_FDR': qcheck.loc[gene, 'FDR'],
            'direction_matches': np.sign(fit.logFC) == expected})
    pd.DataFrame(checks).to_csv(OUT / 'paper_direction_checks.csv', index=False)
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matrix = corr.pivot(index='cohort', columns='baseline_cohort', values='spearman').loc[COHORTS]
    fig, ax = plt.subplots(figsize=(9.5, 3.6), layout='constrained')
    im = ax.imshow(matrix, cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set(xticks=range(8), xticklabels=matrix.columns, yticks=[0, 1],
           yticklabels=['Lrh1 P1 [Nanos3-Cre]', 'Lrh1 P3 [Dhh-Cre]'],
           title='Early Lrh1 perturbations vs frozen eight-cohort effects')
    ax.tick_params(axis='x', rotation=55)
    for i in range(2):
        for j in range(8):
            ax.text(j, i, f'{matrix.iloc[i,j]:.2f}', ha='center', va='center', fontsize=9)
    fig.colorbar(im, ax=ax, shrink=.85, label='Descriptive Spearman')
    fig.savefig(OUT / 'baseline_correlations.png', dpi=180);plt.close(fig)
    pca = pd.read_csv(OUT / 'PCA.csv')
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), layout='constrained')
    for ax, cohort in zip(axes, COHORTS):
        part = pca[pca.cohort == cohort]
        for (group, batch), rows in part.groupby(['group', 'batch']):
            ax.scatter(rows.PC1, rows.PC2, color='#bd5d4d' if group=='KO' else '#246b87',
                marker='o' if batch=='B1' else '^', label=f'{group} / {batch}', s=55)
            for _, r in rows.iterrows():
                ax.annotate(r.gsm[-2:], (r.PC1,r.PC2), xytext=(4,4), textcoords='offset points', fontsize=8)
        ax.set(title=cohort.replace('Lrh1_', '').replace('_', ' '),
               xlabel=f'PC1 ({part.PC1_variance.iloc[0]:.1%})', ylabel=f'PC2 ({part.PC2_variance.iloc[0]:.1%})')
        ax.margins(x=.12, y=.15)
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle('Sample QC: voom expression PCA, not batch-corrected')
    fig.savefig(OUT / 'sample_QC.png', dpi=180);plt.close(fig)
    snap = ROOT / 'followup_snapshot';snap.mkdir(exist_ok=True)
    selected = ['cohort_summary.csv', 'Lrh1_samples.csv', 'Lrh1_mapping_mass.csv',
        'Ppp2ca_unit_audit.csv', 'batch_sensitivity.csv', 'method_sensitivity.csv', 'leave_one_out.csv', 'PCA.csv',
        'baseline_correlations.csv', 'baseline_ranges.csv', 'paper_direction_checks.csv',
        'baseline_correlations.png', 'sample_QC.png']
    selected += [c + '_design.csv' for c in COHORTS]
    for name in selected:
        shutil.copyfile(OUT / name, snap / name)
    sources = ['scripts/fetch_followup.py', 'scripts/run_followup.py', 'scripts/check_followup_reproduction.py',
        'scripts/discover_expansion.py', 'code/prepare_followup.py', 'code/analyze_lrh1.R',
        'code/summarize_followup.py', 'tests/test_followup.py', 'envs/python.sh', 'envs/rscript.sh', 'envs/network.sh']
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    inventory = [{'path': str(p.relative_to(snap)), 'bytes': p.stat().st_size, 'sha256': digest(p)}
                 for p in sorted(snap.rglob('*')) if p.is_file() and p.name!='manifest.json']
    (snap / 'manifest.json').write_text(json.dumps({
        'scope': '16 early Lrh1 bulk samples; Ppp2ca four FPKM profiles audited only; prior snapshots unchanged.',
        'gene_level_counts_preserved_before_DE': True, 'cell_proportions_validated': False,
        'source_sha256': {name: digest(ROOT/name) for name in sources},
        'input_specification_sha256': {name: digest(ROOT/name) for name in ['data/followup_assets.json', 'data/followup_contract.json']},
        'frozen_baseline_effect_sha256': digest(ROOT/'results/all_cohort_comparison/common_logFC.csv.gz'),
        'files': inventory}, indent=2) + '\n')
    print(pd.read_csv(OUT / 'cohort_summary.csv').to_string(index=False))
    print(matrix.to_string())
    print(pd.DataFrame(checks).to_string(index=False))


if __name__ == '__main__':
    main()
