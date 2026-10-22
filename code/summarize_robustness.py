#!/usr/bin/env python3
"""Summarize sensitivity ranges, explicitly not confidence intervals."""
import os
from pathlib import Path

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/robustness'
    metrics=pd.read_csv(out/'sensitivity_metrics.csv')
    summary=metrics.groupby(['gse','kind']).agg(
        variants=('variant','size'),minimum_baseline_rho=('spearman_with_full','min'),
        maximum_baseline_rho=('spearman_with_full','max'),
        minimum_strong_sign_agreement=('sign_agreement_strong','min'))
    summary.to_csv(out/'sensitivity_summary.csv')
    pairs=pd.read_csv(out/'pairwise_sensitivity.csv')
    pairs['pair']=pairs.apply(lambda row:' / '.join(sorted([row.gse,row.other])),axis=1)
    grouped=pairs.groupby(['pair','kind']).agg(variants=('variant','size'),
        minimum_rho=('variant_rho','min'),maximum_rho=('variant_rho','max'))
    grouped['crosses_zero']=(grouped.minimum_rho<0)&(grouped.maximum_rho>0)
    grouped.to_csv(out/'pair_ranges.csv')
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    selected=[('GSE236388','GSE283904'),('GSE236388','GSE294209'),
              ('GSE283904','GSE42601'),('GSE68284','GSE42601'),('GSE101511','GSE236388')]
    fig,ax=plt.subplots(figsize=(9,4.5),layout='constrained')
    labels=[]
    for i,(a,b) in enumerate(selected):
        key=' / '.join(sorted([a,b])); part=pairs[pairs.pair==key]
        baseline=part.baseline_rho.iloc[0]
        ax.plot([part.variant_rho.min(),part.variant_rho.max()],[i,i],lw=3,color='#246b87')
        ax.plot(baseline,i,'o',color='#b45536')
        labels.append(f'{a} / {b}')
    ax.axvline(0,color='gray',lw=.8)
    ax.set(yticks=range(len(labels)),yticklabels=labels,xlim=(-1,1),xlabel='Spearman effect correlation',
           title='One-study-at-a-time sensitivity ranges (not confidence intervals)')
    fig.savefig(out/'key_pair_sensitivity.png',dpi=180)
    plt.close(fig)
    print(summary.to_string())


if __name__=='__main__':
    main()
