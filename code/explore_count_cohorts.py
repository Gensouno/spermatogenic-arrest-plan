#!/usr/bin/env python3
"""Exploratory cross-cohort effects and reference sensitivity; no causal shares."""
import itertools
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from calibrate_reference_batches import design

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--extended-counts', action='store_true')
    args = parser.parse_args()
    out = ROOT / ('results/extended_count_exploration' if args.extended_counts else 'results/count_cohort_exploration')
    out.mkdir(parents=True, exist_ok=True)
    cohorts = pd.read_csv(ROOT / 'results/import_qc/count_cohorts.csv').gse.tolist()
    if args.extended_counts:
        extended = pd.read_csv(ROOT / 'results/import_qc/extended_rnaseq_cohorts.csv')
        extra = extended[extended.unit == 'counts']
        if (extra.min_retained_abundance_fraction < .95).any():
            raise ValueError('Poor count mapping coverage')
        cohorts += extra.gse.tolist()
    effects = {c:pd.read_csv(ROOT / f'results/de/{c}_de.csv').set_index('gene').logFC for c in cohorts}
    # Join only per-cohort effect estimates. Never pool cross-study counts.
    lfc = pd.DataFrame(effects).dropna()
    lfc.to_csv(out / 'shared_logFC.csv.gz')
    spearman = lfc.corr(method='spearman')
    spearman.to_csv(out / 'spearman.csv')
    lfc.corr().to_csv(out / 'pearson.csv')
    references = {name:pd.read_csv(ROOT / f'results/reference_calibration/{name}_cpm.csv.gz',index_col=0).T
                  for name in ['pooled', 'ST_only']}
    rows = []
    for cohort in cohorts:
        x = pd.read_csv(ROOT / f'data/bulk/{cohort}_matrix.tsv.gz',sep='\t',index_col=0)
        meta = pd.read_csv(ROOT / f'data/bulk/{cohort}_samples.tsv', sep='\t').set_index('matrix_sample')
        cpm = x.div(x.sum(axis=0), axis=1)*1e6
        for refname, ref in references.items():
            common = ref.index.intersection(x.index)
            s = ref.loc[common].to_numpy()
            y = cpm.loc[common].to_numpy()
            rng = np.random.default_rng(42)
            order = rng.permutation(len(common))
            tr, te = order[:int(.8*len(common))], order[int(.8*len(common)):]
            for method in ['linear_all','linear_weighted','markers_weighted']:
                use, w = design(s, method)
                fit = np.isin(use,tr)
                use, w = use[fit], w[fit]
                for j, sample in enumerate(x.columns):
                    coef, _ = nnls(s[use]*w[:,None], y[use,j]*w)
                    if coef.sum() <= 0:
                        raise ValueError('Zero coefficients')
                    p = coef/coef.sum()
                    pred = s[te] @ coef
                    obs = y[te,j]
                    raw_r2 = 1-np.sum((obs-pred)**2)/np.sum((obs-obs.mean())**2)
                    logged = np.log2(obs+1)
                    log_r2 = 1-np.sum((logged-np.log2(pred+1))**2)/np.sum((logged-logged.mean())**2)
                    row = {'gse':cohort,'sample':sample,'gsm':meta.loc[sample,'sample_id'],
                           'group':meta.loc[sample,'group'], 'reference':refname,'method':method,
                           'shared_genes':len(common),'fit_genes':len(use),
                           'heldout_cpm_r2':raw_r2,'heldout_log_r2':log_r2}
                    row.update(dict(zip(ref.columns,p)))
                    rows.append(row)
    weights = pd.DataFrame(rows)
    weights.to_csv(out / 'reference_sensitivity.csv',index=False)
    cols = list(next(iter(references.values())).columns)
    group = weights.groupby(['gse','reference','method','group'])[cols].mean()
    delta = group.xs('KO',level='group') - group.xs('WT',level='group')
    delta.to_csv(out / 'KO_minus_WT_weights.csv')
    # Pooled vs ST-only sensitivity for the weighted model; not confidence bounds.
    ranges = delta.xs('linear_weighted',level='method').groupby('gse').agg(['min','max'])
    ranges.to_csv(out / 'weighted_reference_ranges.csv')
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT / '.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(11,4),layout='constrained')
    image = axes[0].imshow(spearman, vmin=-1,vmax=1,cmap='RdBu_r')
    axes[0].set(xticks=range(len(cohorts)),yticks=range(len(cohorts)),
                xticklabels=cohorts,yticklabels=cohorts,title=f'LogFC Spearman; {len(lfc):,} shared tested genes')
    axes[0].tick_params(axis='x',rotation=60)
    for i,j in itertools.product(range(len(cohorts)),repeat=2):
        axes[0].text(j,i,f'{spearman.iloc[i,j]:.2f}',ha='center',va='center',fontsize=9)
    fig.colorbar(image,ax=axes[0],shrink=.7)
    selected = delta.xs('linear_weighted',level='method')['Somatic'].unstack('reference')
    selected.plot.bar(ax=axes[1],color=['#246b87','#b45536'])
    axes[1].axhline(0,color='black',lw=.7)
    axes[1].set(ylabel='KO minus WT fitted RNA weight', title='Somatic weight: reference sensitivity',xlabel=None)
    axes[1].tick_params(axis='x',rotation=60)
    fig.savefig(out / 'exploratory_overview.png',dpi=180)
    plt.close(fig)
    (out / 'interpretation.json').write_text(json.dumps({
        'shared_tested_genes':len(lfc), 'datasets':cohorts, 'genetic_models':len(cohorts)-1,
        'notes':['Rfx2 cohorts are two ages from the same study.',
                 'Sample IDs and genotype evidence audited in results/import_qc/official_sample_evidence.csv.',
                 'RNA weights are model-dependent and not independently validated cell fractions.',
                 'No null-model causal percentage or residual mechanistic DE is inferred.']},indent=2)+'\n')
    print(spearman.to_string())
    print(selected.to_string())


if __name__ == '__main__':
    main()
