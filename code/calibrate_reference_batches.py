#!/usr/bin/env python3
"""Leave one atlas source batch out; evaluate RNA weights, not cell counts."""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[1]
CLASSES = ['Meiotic', 'Postmeiotic', 'Somatic', 'Spermatogonia']


def signature(raw, meta, mask):
    sums = np.stack([raw[mask & (meta.cell_class == c).to_numpy()].sum(axis=0) for c in CLASSES], axis=1)
    if np.any(sums.sum(axis=0) <= 0):
        raise ValueError('Missing reference class')
    return sums / sums.sum(axis=0) * 1e6


def design(s, method):
    use = np.arange(len(s))
    if method == 'markers_weighted':
        # Selection uses training reference only; union of 200 specific genes
        # per class with >=10 CPM. No bulk group labels or test profiles used.
        specificity = s / (s.sum(axis=1, keepdims=True) + 1e-12)
        selected = []
        for k in range(s.shape[1]):
            eligible = np.where(s[:,k] >= 10)[0]
            selected.extend(eligible[np.argsort(specificity[eligible,k])[-200:]])
        use = np.unique(selected)
    weight = np.ones(len(use)) if method == 'linear_all' else 1 / np.sqrt(s[use].mean(axis=1) + 1)
    return use, weight


def estimate(s, y, use, weight):
    coef, _ = nnls(s[use] * weight[:,None], y[use] * weight)
    if coef.sum() <= 0:
        raise ValueError('Zero coefficient total')
    return coef / coef.sum()


def main():
    out = ROOT / 'results/reference_calibration'
    out.mkdir(parents=True, exist_ok=True)
    base = ROOT / 'data/ref/batches'
    data = np.load(base / 'counts.npz')
    raw, genes = data['sums'], data['genes']
    meta = pd.read_csv(base / 'batch_classes.csv')
    eligible = meta.pivot(index='source_batch', columns='cell_class', values='cells').reindex(columns=CLASSES).fillna(0)
    eligible = eligible.index[(eligible >= 20).all(axis=1)].tolist()
    rng = np.random.default_rng(20261022)
    weights = rng.dirichlet(np.ones(4), size=100)
    rows = []
    methods = ['linear_all', 'linear_weighted', 'markers_weighted']
    for batch in eligible:
        heldout = (meta.source_batch == batch).to_numpy()
        train_sig = signature(raw, meta, ~heldout)
        test_sig = signature(raw, meta, heldout)
        for method in methods:
            use, scale = design(train_sig, method)
            for rep, truth in enumerate(weights):
                y = test_sig @ truth
                p = estimate(train_sig, y, use, scale)
                row = {'source_batch':batch, 'method':method, 'replicate':rep,
                       'weight_mae':float(np.abs(p-truth).mean()),
                       'max_weight_error':float(np.abs(p-truth).max()), 'fit_genes':len(use)}
                row.update({f'{c}_error':float(p[k]-truth[k]) for k,c in enumerate(CLASSES)})
                rows.append(row)
        print(f'held out {batch}', flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / 'heldout_mixtures.csv', index=False)
    per_batch = frame.groupby(['method', 'source_batch']).agg(mean_mae=('weight_mae','mean'),
                    median_mae=('weight_mae','median'), max_error=('max_weight_error','max'))
    per_batch.to_csv(out / 'per_batch.csv')
    summary = per_batch.groupby('method').agg(source_batches=('mean_mae','size'),
                    mean_batch_mae=('mean_mae','mean'), worst_batch_mae=('mean_mae','max'))
    summary.to_csv(out / 'summary.csv')
    # Preserve pooled and unenriched-source-only references for real sensitivity.
    for name, mask in [('pooled', np.ones(len(meta),dtype=bool)),
                       ('ST_only', meta.source_batch.str.match(r'^ST[1-6]$').to_numpy())]:
        s = signature(raw, meta, mask)
        pd.DataFrame(s.T, index=CLASSES, columns=genes).to_csv(out / f'{name}_cpm.csv.gz')
    notes = {'seed':20261022, 'min_cells_per_class':20, 'eligible_sources':eligible,
             'mixtures_per_source':100, 'methods':methods,
             'interpretation':'Held-out source-batch generalization within one atlas. Unknown donor independence; no cross-study or cell-number validation. Method comparison is exploratory, not an independent final test.'}
    (out / 'design.json').write_text(json.dumps(notes, indent=2)+'\n')
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10,4), layout='constrained')
    wide = per_batch['mean_mae'].unstack('method').reindex(columns=methods)
    wide.plot.bar(ax=ax, color=['#b45536','#246b87','#6b8e52'])
    ax.set(ylabel='Mean absolute RNA-weight error', xlabel='Held-out source batch',
           title='Cross-source calibration within GSE112393 (not independent animals)')
    ax.legend(title=None, frameon=False)
    fig.savefig(out / 'cross_source_error.png', dpi=180)
    plt.close(fig)
    print(summary.to_string())


if __name__ == '__main__':
    main()
