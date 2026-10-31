#!/usr/bin/env python3
"""Plot all competitive results, keeping model policies and coverage explicit."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from run_independent_contexts import ROOT,COHORTS


def plot(out):
    d=pd.read_csv(out/'camera_results.csv');fig,axes=plt.subplots(1,2,figsize=(12,5.1))
    colors={'estimated':'#174c74','fixed_0p01':'#b85a20'}
    labels={'estimated':'Estimated correlation (main)','fixed_0p01':'Assumed correlation 0.01 (sensitivity)'}
    for policy,offset in [('estimated',-.1),('fixed_0p01',.1)]:
        z=d[d.policy==policy].set_index('cohort').loc[COHORTS]
        for ax,values in zip(axes,[z.correlation_raw_or_assumed,-np.log10(z.BY)]):
            ax.scatter(values,np.arange(4)+offset,color=colors[policy],label=labels[policy],s=45,zorder=3)
    axes[0].set_yticks(range(4),['Rif1 P14 (30/36)','Hormad SFB P19 (36/36)','Mael P16 (23/36)','Morc1 P14.5 (31/36)'])
    axes[1].set_yticks(range(4),[])
    for ax in axes:
        ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
    axes[0].set_xlabel('Estimated or assumed residual correlation')
    axes[1].set_xlabel('-log10(BY-adjusted two-sided P value)')
    axes[1].axvline(-np.log10(.05),color='#555555',ls='--',lw=1)
    axes[1].text(-np.log10(.05)+.12,.98,'BY = 0.05',fontsize=8,
        transform=axes[1].get_xaxis_transform(),va='top')
    axes[0].set_title('Correlation policy and measured priorities',fontsize=10)
    axes[1].set_title('Eight tests adjusted together; all directions Down',fontsize=10)
    handles,legend_labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,legend_labels,loc='lower center',bbox_to_anchor=(.5,.09),ncol=2,fontsize=9)
    fig.text(.02,.03,'Estimated-policy df = 4; fixed-policy df = universe size - 2. Policies differ in both variance inflation and test df.\nPost-hoc set association; no individual-gene discoveries or composition adjustment. Original primary family unchanged.',fontsize=8)
    fig.tight_layout(rect=(0,.18,1,.98));fig.savefig(out/'correlation_sensitivity.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results/competitive_set')
    plot(parser.parse_args().output)
