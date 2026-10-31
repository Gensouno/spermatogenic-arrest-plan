#!/usr/bin/env python3
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from run_stage_specificity import ROOT,COHORTS


def main(out):
    data=pd.read_csv(out/'specificity_summary.csv');fig,axes=plt.subplots(2,2,figsize=(11,8))
    colors={'priority':'#222222','abundance_only':'#e69f00','stage_abundance_dominance':'#0072b2'}
    for ax,cohort in zip(axes.flat,COHORTS):
        ticklabels=[]
        for j,ref in enumerate(['Chen','Green']):
            d=data[(data.cohort==cohort)&(data.reference==ref)&(data.abundance_caliper==1)]
            n=int(d.matched_priorities.iloc[0]);target=d.priority_mean_logFC.iloc[0]
            ticklabels.append(f'{ref}\n{n}/36 matched')
            ax.plot(target,j,'D',color=colors['priority'],label='Fixed priorities' if j==0 else None)
            for scheme,dy,marker,label in [('abundance_only',-.12,'o','Abundance-matched background'),('stage_abundance_dominance',.12,'s','Stage/abundance/dominance background')]:
                value=d[d.comparison==scheme].background_mean_logFC.iloc[0]
                ax.plot([value,target],[j+dy,j],color=colors[scheme],alpha=.3)
                ax.plot(value,j+dy,marker,color=colors[scheme],label=label if j==0 else None)
        ax.axvline(0,color='gray',lw=.7);ax.set(yticks=[0,1],yticklabels=ticklabels,title=cohort,
            xlabel='Equal-priority mean log2(perturbation/control)',ylim=(-.5,1.8));ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=1,fontsize=9,bbox_to_anchor=(.5,.02))
    fig.suptitle('Reference-informed descriptive matching at the primary abundance caliper',fontsize=12)
    fig.text(.02,.005,'Same eligible targets within each view; reused background genes are not independent animals. No P values or confidence intervals.',fontsize=8)
    fig.tight_layout(rect=(0,.16,1,.95));fig.savefig(out/'matched_backgrounds.png',dpi=180);plt.close(fig)
    resolution=pd.read_csv(out/'resolution_summary.csv');fig,axes=plt.subplots(2,2,figsize=(11,8))
    for ax,cohort in zip(axes.flat,COHORTS):
        ticks=[]
        for j,ref in enumerate(['Chen','Green']):
            d=resolution[(resolution.cohort==cohort)&(resolution.reference==ref)&(resolution.abundance_caliper==1)]
            ticks.append(f'{ref}\n{int(d.matched_priorities.iloc[0])}/36 matched')
            for k,(scheme,color,marker,label) in enumerate([('abundance_only','#e69f00','o','Abundance only'),('stage_abundance_dominance','#0072b2','s','Add group-mean peak + dominance'),('stage_peak_abundance_dominance','#cc79a7','^','Also require single-stage peak')]):
                value=d[d.comparison==scheme].mean_logFC_difference.iloc[0]
                ax.plot(value,j+(k-1)*.13,marker,color=color,label=label if j==0 else None)
        ax.axvline(0,color='gray',lw=.7);ax.set(yticks=[0,1],yticklabels=ticks,title=cohort,
            xlabel='Mean priority minus matched-background log2 effect',ylim=(-.5,1.5));ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='lower center',ncol=1,fontsize=9,bbox_to_anchor=(.5,.02))
    fig.suptitle('Additional peak-stage constraint: compare schemes on common targets',fontsize=12)
    fig.text(.02,.005,'Different expression descriptors need not be nested. These are post-hoc descriptive comparisons, without P values or confidence intervals.',fontsize=8)
    fig.tight_layout(rect=(0,.16,1,.95));fig.savefig(out/'stage_resolution.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default=str(ROOT/'results/stage_specificity'));a=p.parse_args();main(Path(a.output))
