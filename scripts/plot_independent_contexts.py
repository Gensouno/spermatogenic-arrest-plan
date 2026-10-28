#!/usr/bin/env python3
"""Static scientific QC and fixed-priority context figures."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from run_independent_contexts import ROOT,COHORTS

p=argparse.ArgumentParser();p.add_argument('--output',default=str(ROOT/'results/independent_contexts'));args=p.parse_args()
out=Path(args.output)
d=pd.read_csv(out/'priority_stability.csv');genes=sorted(d.symbol.unique())
values=d.pivot(index='symbol',columns='cohort',values='full_logFC').reindex(index=genes,columns=COHORTS)
fig,ax=plt.subplots(figsize=(9,11.5));cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd')
im=ax.imshow(values.to_numpy(),aspect='auto',cmap=cmap,vmin=-4,vmax=4)
ax.set_xticks(range(4),['Rif1 P14\nphenotype unresolved','Hormad1/2 SFB P19','Mael P16\nHET controls','Morc1 P14.5\nHET controls'],fontsize=9)
ax.set_yticks(range(len(genes)),genes,fontsize=8)
for r,g in enumerate(genes):
    for c,cohort in enumerate(COHORTS):
        row=d[(d.symbol==g)&(d.cohort==cohort)].iloc[0]
        if row.all_down_1p5fold:ax.text(c,r,'*',ha='center',va='center',color='white' if abs(row.full_logFC)>2.5 else 'black',fontsize=12)
fig.colorbar(im,ax=ax,label='Full-model log2(perturbation/control); color clipped at +/-4',shrink=.65)
ax.set_title('Frozen descriptive priorities in separate source contexts',fontsize=11)
fig.text(.02,.02,'* >=1.5-fold down in full + all six single-library omissions; descriptive only.\nGray: filtered, absent or nonunique. No added primary votes; animal independence unverified.',fontsize=8)
fig.tight_layout(rect=(0,.06,1,1));fig.savefig(out/'priority_context.png',dpi=180);plt.close(fig)
p=pd.read_csv(out/'PCA.csv');fig,axes=plt.subplots(1,4,figsize=(17,4.5))
for ax,cohort in zip(axes,COHORTS):
    s=p[p.cohort==cohort]
    for group,color in [('Control','#377eb8'),('Perturbation','#e41a1c')]:
        t=s[s.group==group];ax.scatter(t.PC1,t.PC2,label=group,color=color)
        for row in t.itertuples():
            offset=(4,-14) if cohort=='Hormad_dsfb_P19' and row.matrix_sample=='dsfb_1' else (4,8) if cohort=='Hormad_dsfb_P19' and row.matrix_sample=='dsfb_2' else (3,3)
            ax.annotate(row.matrix_sample.replace('geneCount_mrna_',''),(row.PC1,row.PC2),fontsize=7,xytext=offset,textcoords='offset points')
    ax.margins(x=.18,y=.18)
    ax.set(title=cohort,xlabel=f'PC1 ({s.PC1_variance.iloc[0]:.1%})',ylabel=f'PC2 ({s.PC2_variance.iloc[0]:.1%})')
axes[0].legend(fontsize=8);fig.suptitle('Voom-expression PCA; no sample excluded using this diagnostic')
fig.tight_layout();fig.savefig(out/'sample_PCA.png',dpi=180);plt.close(fig)
