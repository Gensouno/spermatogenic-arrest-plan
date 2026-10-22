#!/usr/bin/env python3
"""Compare log-space fitting with linear RNA-mixture fitting; no cell-count claims."""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'mixture_audit'


def fit_predict(s, y, train, mode):
    a = np.log2(s + 1) if mode == 'legacy_log' else s
    b = np.log2(y + 1) if mode == 'legacy_log' else y
    coef, _ = nnls(a[train], b[train])
    total = coef.sum()
    if total <= 0:
        raise ValueError('Zero mixture coefficient sum')
    p = coef / total
    # Legacy script normalizes coefficients before predicting. Linear NNLS
    # retains its fitted scale for held-out reconstruction; weights sum to one.
    pred = a @ (p if mode == 'legacy_log' else coef)
    return p, pred, b, total


def r2(obs, pred):
    denom = np.sum((obs - obs.mean()) ** 2)
    return float(1 - np.sum((obs - pred) ** 2) / denom) if denom > 0 else None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Regenerated from frozen raw atlas assets by reference_calibrate; do not
    # depend on the old untracked, rounded signature file.
    sig = ROOT / 'results/reference_calibration/pooled_cpm.csv.gz'
    bulk = ROOT / 'data/raw/GSE236388/GSE236388_XAP5_counts.txt.gz'
    s = pd.read_csv(sig, index_col=0).T
    x = pd.read_csv(bulk, sep='\t', index_col=0).groupby(level=0).sum()
    x = x.div(x.sum(axis=0), axis=1) * 1e6
    genes = s.index.intersection(x.index)
    a = s.loc[genes].to_numpy(dtype=float)
    y = x.loc[genes].to_numpy(dtype=float)
    if not (np.isfinite(a).all() and np.isfinite(y).all()):
        raise ValueError('Non-finite input')
    if a.min() < 0 or y.min() < 0:
        raise ValueError('Negative input')
    rng = np.random.default_rng(42)
    split = rng.permutation(len(genes))
    train = split[:int(.8 * len(genes))]
    test = split[int(.8 * len(genes)):]
    pd.DataFrame({'gene': genes, 'partition': np.where(np.isin(np.arange(len(genes)), train), 'fit', 'heldout')}).to_csv(OUT / 'gene_split.csv', index=False)
    simulations = []
    # Exact-reference mixtures are an algebra sanity check, not independent
    # biological validation. Compare both models on exactly the same mixtures.
    for rep, p in enumerate(rng.dirichlet(np.ones(a.shape[1]), size=200)):
        synthetic = a @ p
        for mode in ('legacy_log', 'linear_cpm'):
            est, pred, obs, scale = fit_predict(a, synthetic, train, mode)
            simulations.append({'replicate': rep, 'method': mode,
                                'mean_absolute_weight_error': np.abs(est-p).mean(),
                                'max_absolute_weight_error': np.abs(est-p).max(),
                                'heldout_r2': r2(obs[test], pred[test])})
    sim = pd.DataFrame(simulations)
    sim.to_csv(OUT / 'exact_reference_calibration.csv', index=False)
    rows = []
    for j, name in enumerate(x.columns):
        for mode in ('legacy_log', 'linear_cpm'):
            p, pred, obs, scale = fit_predict(a, y[:, j], train, mode)
            row = dict(sample=name, group='KO' if 'cKO' in name else 'WT', method=mode,
                       heldout_r2=r2(obs[test], pred[test]), coefficient_sum=scale)
            row.update(dict(zip(s.columns, p)))
            # Same scoring scale is necessary for comparisons across methods.
            # Convert linear predictions to log only for this additional score.
            log_pred = pred if mode == 'legacy_log' else np.log2(pred + 1)
            row['heldout_log_r2'] = r2(np.log2(y[test, j] + 1), log_pred[test])
            rows.append(row)
    real = pd.DataFrame(rows)
    real.to_csv(OUT / 'real_heldout.csv', index=False)
    delta = real.groupby(['method', 'group'])[list(s.columns)].mean()
    delta.to_csv(OUT / 'group_weights.csv')
    summary = sim.groupby('method').agg(
        median_weight_mae=('mean_absolute_weight_error', 'median'),
        maximum_weight_error=('max_absolute_weight_error', 'max'),
        median_heldout_r2=('heldout_r2', 'median'))
    summary.to_csv(OUT / 'calibration_summary.csv')
    assert sim.loc[sim.method == 'linear_cpm', 'max_absolute_weight_error'].max() < 1e-6
    provenance = {'seed': 42, 'shared_genes': len(genes), 'fit_genes': len(train),
                  'heldout_genes': len(test), 'simulated_mixtures': 200,
                  'scope': 'Exact-reference algebra calibration and real-data reconstruction; not cell-fraction validation',
                  'inputs': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (sig, bulk)}}
    (OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    os.environ.setdefault('MPLCONFIGDIR', str(OUT / '.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    methods = ['legacy_log', 'linear_cpm']
    axes[0].boxplot([sim.loc[sim.method == m, 'mean_absolute_weight_error'] for m in methods],
                    tick_labels=['Legacy log', 'Linear CPM'], showfliers=False)
    axes[0].set(ylabel='Mean absolute RNA-weight error', title='Known mixtures (same reference)')
    for m, color in zip(methods, ['#b45536', '#246b87']):
        part = real[real.method == m]
        axes[1].plot(range(6), part.heldout_log_r2, 'o-', label=m, color=color)
    axes[1].set(xticks=range(6), xticklabels=['WT1', 'WT2', 'WT3', 'KO1', 'KO2', 'KO3'],
                ylabel='Held-out R² on log2(CPM+1)', title='Real samples: common scoring scale')
    axes[1].legend(frameon=False)
    fig.savefig(OUT / 'mixture_audit.png', dpi=180)
    plt.close(fig)
    print(summary.to_string())
    print(real[['sample', 'method', 'heldout_r2', 'heldout_log_r2']].to_string(index=False))


if __name__ == '__main__':
    main()
