#!/usr/bin/env python3
"""Controlled broad-vs-fine atlas check, not cell proportion validation.

The 12 published germ centroids are averages of per-cell normalized UMI,
not count pseudobulks. Both resolutions use these SAME centroids and the
same somatic profile, genes, weights and gene folds. External labels are
not used to select genes, tune parameters or relabel the atlas clusters.
"""
import hashlib
import json
import os
from pathlib import Path
import urllib.request

import numpy as np
import pandas as pd
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[1]
FILENAME = 'GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx'
URL = f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112393/suppl/{FILENAME}'
PAPER = 'https://pubmed.ncbi.nlm.nih.gov/30146481/'
GROUPS = {'Spermatogonia': ['GC1'],
          'Meiotic': [f'GC{i}' for i in range(2, 9)],
          'Postmeiotic': [f'GC{i}' for i in range(9, 13)]}


def read_centroids(path):
    average = pd.read_excel(path, sheet_name='Average', header=5, index_col=0)
    total = pd.read_excel(path, sheet_name='Total', header=4, index_col=0)
    columns = [f'GC{i}' for i in range(1, 13)]
    assert list(average.columns) == columns == list(total.columns)
    cells = average.loc['nCell'].astype(int)
    assert cells.sum() == 20646
    assert np.array_equal(cells, total.loc['cell'].astype(int))
    average = average.drop(index=['nGene', 'nCell']).astype(float)
    total = total.drop(index=['nGene', 'cell']).astype(float)
    assert average.index.equals(total.index) and average.shape == (24475, 12)
    assert np.isfinite(average.to_numpy()).all() and (average >= 0).all().all()
    # Verify the logarithm base from TWO published sheets, not convention.
    # Total = ln(1 + nCell * expm1(Average)) for the natural-log scale.
    error = float(np.abs(np.log1p(np.expm1(average) * cells) - total).to_numpy().max())
    if error > 1e-8:
        raise ValueError(f'Average/Total natural-log identity failed: {error}')
    linear = np.expm1(average)
    if not np.allclose(linear.sum(), 10000, atol=1e-6):
        raise ValueError('Expected each mean normalized profile to sum to 10,000')
    invalid = np.array([not isinstance(gene, str) for gene in linear.index])
    excluded = linear.loc[invalid].copy()
    excluded.index = excluded.index.map(str)
    linear = linear.loc[~invalid].copy()
    if linear.index.has_duplicates:
        raise ValueError('Duplicate valid gene symbols require explicit resolution')
    return linear * 100, cells, error, excluded


def collapse_centroids(fine, cells):
    broad = pd.DataFrame(index=fine.index)
    for group, clusters in GROUPS.items():
        broad[group] = fine[clusters].mul(cells[clusters]).sum(axis=1) / cells[clusters].sum()
    return broad


def fit(reference, y, scale):
    if (reference.ndim != 2 or y.shape != (len(reference),) or scale.shape != y.shape
            or not np.isfinite(reference).all() or not np.isfinite(y).all()
            or not np.isfinite(scale).all() or np.any(reference < 0)
            or np.any(y < 0) or np.any(scale <= 0)):
        raise ValueError('Invalid nonnegative abundance or gene weights')
    coef, _ = nnls(reference * scale[:, None], y * scale, maxiter=1000)
    if coef.sum() <= 0:
        raise ValueError('Empty nonnegative fit')
    return coef


def main():
    path = ROOT / 'data/raw/ref' / FILENAME
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(URL, timeout=45) as response:
            data = response.read()
        if not data.startswith(b'PK'):
            raise ValueError('Atlas supplement is not an XLSX archive')
        path.write_bytes(data)
    fine, cells, error, excluded = read_centroids(path)
    broad = collapse_centroids(fine, cells)
    external = pd.read_csv(ROOT / 'results/external_sorted_reference/source_RPKM.csv', index_col=0)
    pooled = pd.read_csv(ROOT / 'results/reference_calibration/pooled_cpm.csv.gz', index_col=0).T
    genes = fine.index.intersection(pooled.index).intersection(external.index).sort_values()
    assert len(genes) > 10000
    # The selected-gene universe is renormalized identically in every model.
    fine = fine.loc[genes]
    broad = broad.loc[genes]
    external = external.loc[genes]
    y = external / external.sum() * 1e6
    out = ROOT / 'results/stage_reference'
    out.mkdir(parents=True, exist_ok=True)
    excluded.to_csv(out / 'excluded_excel_date_gene_rows.csv')
    fine.to_csv(out / 'germ_centroids_linear_cpm.csv.gz')
    cells.rename('cells').to_csv(out / 'germ_cluster_sizes.csv')
    pd.DataFrame({'gene': genes, 'fold': np.random.default_rng(112393).permutation(len(genes)) % 5}).to_csv(out / 'gene_folds.csv', index=False)
    folds = pd.read_csv(out / 'gene_folds.csv')['fold'].to_numpy()
    full_weights, fold_weights, metrics, predictions = [], [], [], []
    for somatic_source in ['pooled', 'ST_only']:
        original = pd.read_csv(ROOT / f'results/reference_calibration/{somatic_source}_cpm.csv.gz', index_col=0).T.loc[genes]
        detailed = fine.copy()
        detailed.insert(0, 'Somatic', original['Somatic'])
        collapsed = broad.copy()
        collapsed.insert(0, 'Somatic', original['Somatic'])
        models = {'original_4': original, 'centroid_4': collapsed, 'centroid_13': detailed}
        models = {name: ref / ref.sum() * 1e6 for name, ref in models.items()}
        # Keep a single gene weighting metric for all three model resolutions.
        scale = 1 / np.sqrt(models['centroid_4'].mean(axis=1).to_numpy() + 1)
        for name, reference in models.items():
            reference.to_csv(out / f'{somatic_source}_{name}_reference.csv.gz')
            s = reference.to_numpy()
            condition = float(np.linalg.cond(s * scale[:, None]))
            for stage in external.columns:
                target = y[stage].to_numpy()
                coef = fit(s, target, scale)
                p = coef / coef.sum()
                aggregate = {}
                for group, clusters in GROUPS.items():
                    aggregate[group] = float(sum(p[list(reference.columns).index(c)] for c in clusters)) if name == 'centroid_13' else float(p[list(reference.columns).index(group)])
                row = {'somatic_source': somatic_source, 'model': name, 'stage': stage,
                       'shared_genes': len(genes), 'largest_component': reference.columns[np.argmax(p)],
                       'coefficient_sum': float(coef.sum()), 'weighted_condition_number': condition,
                       'Somatic': float(p[list(reference.columns).index('Somatic')]), **aggregate}
                if name == 'centroid_13':
                    row.update(dict(zip(reference.columns, p)))
                    row['GC2_GC3_transition'] = float(p[2:4].sum())
                full_weights.append(row)
                prediction = np.empty_like(target)
                for fold in range(5):
                    train, test = folds != fold, folds == fold
                    assert not np.any(train & test) and np.all(train | test)
                    fold_coef = fit(s[train], target[train], scale[train])
                    prediction[test] = s[test] @ fold_coef
                    for component, value in zip(reference.columns, fold_coef / fold_coef.sum()):
                        fold_weights.append({'somatic_source': somatic_source, 'model': name,
                            'stage': stage, 'fold': fold, 'component': component, 'RNA_weight': float(value)})
                metrics.append({'somatic_source': somatic_source, 'model': name, 'stage': stage,
                    'shared_genes': len(genes), 'crossfit_log2_RMSE': float(np.sqrt(np.mean((np.log2(target + .5) - np.log2(prediction + .5)) ** 2))),
                    'crossfit_weighted_RMSE': float(np.sqrt(np.mean(((target - prediction) * scale) ** 2)))})
                predictions.append(pd.DataFrame({'gene': genes, 'fold': folds,
                    'somatic_source': somatic_source, 'model': name, 'stage': stage,
                    'observed_relative_RPKM': target, 'crossfit_predicted_abundance': prediction}))
    weights = pd.DataFrame(full_weights)
    weights.to_csv(out / 'fitted_weights.csv', index=False)
    pd.DataFrame(fold_weights).to_csv(out / 'fold_weights.csv', index=False)
    pd.DataFrame(metrics).to_csv(out / 'crossfit_metrics.csv', index=False)
    pd.concat(predictions).to_csv(out / 'crossfit_predictions.csv.gz', index=False)
    audit = {'source_url': URL, 'paper': PAPER, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'natural_log_Average_Total_identity_max_error': error,
        'normalization': 'exp(Average)-1 sums to 10,000 in each cluster; multiply by 100 before shared-gene renormalization.',
        'germ_cells': int(cells.sum()), 'original_genes': 24475,
        'excluded_nonstring_gene_rows': len(excluded), 'shared_genes': len(genes),
        'cluster_annotation': {'GC1': 'Spermatogonia', 'GC2-GC3': 'Rare transition into meiosis; not asserted exact LZ identity',
            'GC4-GC8': 'Spermatocytes; not assigned exact prophase substages',
            'GC9-GC11': 'Round spermatids', 'GC12': 'Elongated spermatids'},
        'design': 'Compare same-table cell-number-weighted collapsed centroids vs 12 germ centroids plus one somatic profile; same genes, weights and five gene folds. Original pseudobulk 4-class reference is an additional comparator.',
        'limitations': ['Exploratory refinement motivated by failure on these SAME external profiles; not an untouched confirmatory validation set.',
            'Averages of per-cell normalized counts are not sum-count pseudobulks. Somatic profile comes from a different aggregation scheme of the same atlas.',
            'Early transition clusters contain only 132 and 61 cells across the atlas; per-batch fine cluster labels unavailable here.',
            'Gene folds are not independent biological samples; correlated genes remain. No claim of donor or sample holdout.',
            'Seven malformed external rows remain excluded. RPKM vs UMI, ages and selected-gene universe remain confounded.',
            'Atlas gene symbols converted to Excel dates are recorded and excluded, not guessed or silently merged.',
            'No full known-proportion mixture, absolute cell counts or causal composition percentages are validated.']}
    (out / 'provenance.json').write_text(json.dumps(audit, indent=2) + '\n')
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    for ax, model in zip(axes, ['centroid_4', 'centroid_13']):
        plot = weights[(weights.somatic_source == 'pooled') & (weights.model == model)].set_index('stage')
        plot[['Somatic', *GROUPS]].plot.bar(stacked=True, ax=ax, color=['#246b87', '#e0a458', '#8a5a83', '#6b8e52'])
        ax.set(xlabel=None, ylabel='Fitted RNA weights (not cell fractions)', title=f'Same centroid table: {model}')
        ax.tick_params(axis='x', rotation=0)
        ax.get_legend().remove()
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=4, frameon=False)
    fig.savefig(out / 'controlled_resolution_comparison.png', dpi=180)
    plt.close(fig)
    print(weights[['somatic_source', 'model', 'stage', 'largest_component', *GROUPS]].to_string(index=False))
    print(pd.DataFrame(metrics).to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
