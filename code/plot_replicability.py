#!/usr/bin/env python3
"""Descriptive plots with all stability priorities and explicit discovery status."""
import os
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from annotate_replicability import CHEN
from analyze_replicability import PRIMARY


def main():
    out=ROOT/'results/replicability'
    c=pd.read_csv(out/'candidate_stage_context.csv')
    stable=pd.read_csv(out/'stable_descriptive_priorities.csv').sort_values(['Chen_peak_group','gene'])
    formal=pd.read_csv(out/'primary_BY_findings.csv')
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    values=[int((c.descriptive_direction=='up').sum()),int((c.descriptive_direction=='down').sum()),
            int((stable.descriptive_direction=='up').sum()),int((stable.descriptive_direction=='down').sum()),len(formal)]
    ax.bar(['Descriptive\nup','Descriptive\ndown','LOO retained\nup','LOO retained\ndown','Joint BY\nhypotheses'],values,
           color=['#bd6552','#527b96','#bd6552','#527b96','#777777'])
    for i,n in enumerate(values):ax.text(i,n+1,str(n),ha='center')
    ax.set(ylabel='Genes (first 4 bars); hypotheses (last bar)',ylim=(0,max(values)*1.18),
           title='1.5-fold in >=3 of 4 genetic units; no formal discoveries')
    fig.savefig(out/'descriptive_support.png',dpi=180);plt.close(fig)
    transfer=pd.read_csv(out/'external_context_summary.csv')
    transfer=transfer[(transfer.subset=='primary_LOO_priorities') & (transfer.candidate_direction=='all')]
    fig,ax=plt.subplots(figsize=(10,4.5),layout='constrained')
    categories={'Retained >=1.5-fold in all LOO':transfer.baseline_and_all_LOO_1p5fold.to_numpy(),
        'Baseline >=1.5-fold, not retained':(transfer.directional_1p5fold-transfer.baseline_and_all_LOO_1p5fold).to_numpy(),
        'Tested, below or opposite threshold':(transfer.tested_exact_unique_symbol-transfer.directional_1p5fold).to_numpy(),
        'Not tested / unmapped / ambiguous':(transfer.fixed_genes-transfer.tested_exact_unique_symbol).to_numpy()}
    left=np.zeros(len(transfer))
    for (label,values),colour in zip(categories.items(),['#416f8a','#9db7c5','#c7c1ba','#eeeeee']):
        ax.barh(range(len(transfer)),values,left=left,color=colour,label=label,edgecolor='white')
        left+=values
    ax.set_yticks(range(len(transfer)),['Dcaf2 P9','Pramef12 P2 [HET control]',
        'Pramef12 P7 [HET control]','Msl3 P45 [phenotype-negative]',
        'Lrh1 P1 [Nanos3-Cre]','Lrh1 P3 [Dhh-Cre]'],fontsize=8)
    ax.invert_yaxis();ax.set(xlim=(0,36),xlabel='Fixed 36 downward priorities; missing retained in denominator',
        title='Post-hoc contexts: different ages and controls; not untouched validation')
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.17),ncol=2,frameon=False,fontsize=8)
    fig.savefig(out/'external_context_coverage.png',dpi=180);plt.close(fig)
    if stable.empty:return
    profiles=pd.read_csv(ROOT/'results/reference_limits/chen_stage_profiles.csv.gz',index_col='gene')
    stages=[s for group in CHEN.values() for s in group]
    v=np.log2(profiles.reindex(stable.gene)[[f'Chen:{s}:equal_cell_mean' for s in stages]].to_numpy()+.5)
    z=np.full_like(v,np.nan);measured=np.isfinite(v).all(axis=1)
    present=v[measured];sd=present.std(axis=1,keepdims=True)
    z[measured]=np.divide(present-present.mean(axis=1,keepdims=True),sd,
                         out=np.full_like(present,np.nan),where=sd>1e-12)
    effects=stable[[g+'_logFC' for g in PRIMARY]].to_numpy()
    fig,axes=plt.subplots(1,2,figsize=(12,max(7,len(stable)*.23)),gridspec_kw={'width_ratios':[1,3]},layout='constrained')
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd')
    a=axes[0].imshow(effects,aspect='auto',cmap=cmap,vmin=-2,vmax=2)
    b=axes[1].imshow(z,aspect='auto',cmap=cmap,vmin=-2,vmax=2)
    axes[0].set_yticks(range(len(stable)),stable.gene,fontsize=7)
    axes[1].set_yticks([])
    axes[0].set_xticks(range(5),['XAP5','N6amt1','Rfx2\nP29-30','Rfx2\nP21','Hnrnpgt'],rotation=60,ha='right',fontsize=8)
    axes[1].set_xticks(range(len(stages)),stages,rotation=60,ha='right',fontsize=8)
    axes[0].set_title('KO - WT log2FC\ncolour clipped at +/-2',fontsize=10)
    axes[1].set_title('Chen equal-cell CPM: within-gene z of log2(CPM+0.5)\nExpression context; not cell fractions',fontsize=10)
    fig.colorbar(a,ax=axes[0],shrink=.4,label='log2FC');fig.colorbar(b,ax=axes[1],shrink=.4,label='within-gene z')
    fig.suptitle(f'All {len(stable)} descriptive LOO-stable priorities; none are formal discoveries',fontsize=11)
    fig.savefig(out/'stable_priority_context.png',dpi=180);plt.close(fig)


if __name__=='__main__':main()
