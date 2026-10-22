#!/usr/bin/env python3
"""Compare effects across six cohorts; retain expression-unit distinctions."""
import json
import argparse
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--arrays',action='store_true',help='Include microarray effect estimates as sensitivity')
    args=parser.parse_args()
    primary=pd.read_csv(ROOT/'results/import_qc/count_cohorts.csv').gse.tolist()
    extended=pd.read_csv(ROOT/'results/import_qc/extended_rnaseq_cohorts.csv')
    cohorts=primary+extended.gse.tolist()
    array_cohorts=[]
    if args.arrays:
        array_cohorts=pd.read_csv(ROOT/'results/import_qc/microarray_cohorts.csv').gse.tolist()
        cohorts+=array_cohorts
    effects={c:pd.read_csv(ROOT/f'results/de/{c}_de.csv').set_index('gene').logFC for c in cohorts}
    all_effects=pd.DataFrame(effects)
    common=all_effects.dropna()
    out=ROOT/('results/all_cohort_comparison' if args.arrays else 'results/rnaseq_effect_comparison')
    out.mkdir(parents=True,exist_ok=True)
    common.to_csv(out/'common_logFC.csv.gz')
    corr=common.corr(method='spearman')
    corr.to_csv(out/'common_gene_spearman.csv')
    paired=all_effects.corr(method='spearman',min_periods=2000)
    paired.to_csv(out/'pairwise_gene_spearman.csv')
    present=all_effects.notna().astype(int)
    (present.T @ present).to_csv(out/'pairwise_gene_numbers.csv')
    # Markers are a fixed diagnostic list, not a data-selected discovery panel.
    marker_list=['Ddx4','Dazl','Zbtb16','Gfra1','Kit','Stra8','Sycp3','Dmc1','Spo11',
                 'Prm1','Prm2','Tnp1','Tnp2','Sox9','Wt1','Cyp11a1','Star','Rfx2','N6amt1','Hnrnpgt']
    all_effects.reindex(marker_list).to_csv(out/'diagnostic_markers.csv')
    # PCA: descriptive QC in each study's own transformed expression space.
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    nrows=math.ceil(len(cohorts)/3)
    fig,axes=plt.subplots(nrows,3,figsize=(12,3.5*nrows),layout='constrained')
    label_rows=[]
    for ax,c in zip(axes.flat,cohorts):
        pca=pd.read_csv(ROOT/f'results/import_qc/{c}_pca.csv')
        for group,color in [('WT','#246b87'),('KO','#b45536')]:
            part=pca[pca.group==group]
            ax.scatter(part.PC1,part.PC2,label=group,color=color)
            for index,(_,row) in enumerate(part.iterrows(),1):
                label=f'{group}{index}'
                label_rows.append({'gse':c,'plot_label':label,'matrix_sample':row['sample']})
                ax.annotate(label, (row.PC1,row.PC2),fontsize=7,
                            xytext=(3,3),textcoords='offset points')
        ax.margins(.15)
        ax.set(title=c+(' (array)' if c in array_cohorts else ' (TPM)' if c=='GSE294209' else ' (counts)'),
               xlabel=f'PC1 ({pca.PC1_variance.iloc[0]:.0%})',
               ylabel=f'PC2 ({pca.PC2_variance.iloc[0]:.0%})')
        ax.legend(frameon=False,fontsize=8)
    for ax in axes.flat[len(cohorts):]:
        ax.axis('off')
    fig.savefig(out/'within_study_PCA.png',dpi=170)
    plt.close(fig)
    pd.DataFrame(label_rows).to_csv(out/'PCA_sample_labels.csv',index=False)
    fig,ax=plt.subplots(figsize=(7,6),layout='constrained')
    im=ax.imshow(corr,cmap='RdBu_r',vmin=-1,vmax=1)
    labels=[c+(' [array]' if c in array_cohorts else ' [TPM]' if c=='GSE294209' else '') for c in cohorts]
    ax.set(xticks=range(len(cohorts)),yticks=range(len(cohorts)),xticklabels=labels,
           yticklabels=labels,title=f'Spearman logFC: {len(common):,} shared tested genes')
    ax.tick_params(axis='x',rotation=65)
    for i in range(len(cohorts)):
        for j in range(len(cohorts)):
            ax.text(j,i,f'{corr.iloc[i,j]:.2f}',ha='center',va='center',fontsize=9)
    fig.colorbar(im,ax=ax,shrink=.7)
    fig.savefig(out/('eight_cohort_effects.png' if args.arrays else 'six_cohort_effects.png'),dpi=180)
    plt.close(fig)
    notes={'datasets':cohorts,'common_tested_genes':len(common),'genetic_perturbations':len(cohorts)-1,
           'arrays_are_sensitivity_only':array_cohorts,
           'TPM_is_sensitivity_only':'GSE294209',
           'caveats':['GSE101511 annotation release 79 chosen for compatibility; original GTF unspecified.',
                      'TPM logFC is not quantitatively equivalent to TMM-normalized count logFC.',
                      'Correlations are descriptive; gene counts do not create independent biological replicates.',
                      'Common-gene and pairwise-gene correlations use different gene universes.']}
    (out/'design.json').write_text(json.dumps(notes,indent=2)+'\n')
    print(corr.to_string())
    print('Common genes:',len(common))


if __name__=='__main__':
    main()
