#!/usr/bin/env python3
"""Compact diagnostic figures and machine-readable decision evidence."""
import json

import numpy as np
import pandas as pd

from diagnose_reference_limits import ROOT


def main():
    out=ROOT/'results/reference_limits';technical=pd.read_csv(out/'technical_comparison.csv')
    transfer=pd.read_csv(out/'transfer_comparison.csv');held=pd.read_csv(out/'source_library_holdout.csv')
    focal=technical[(technical.model=='centroid_13')&technical.author_type.isin(['Leptotene','Zygotene','L','Z'])]
    focal.to_csv(out/'technical_focal_summary.csv',index=False)
    tf=transfer[transfer.author_type.isin(['Leptotene','Zygotene'])]
    tf.to_csv(out/'transfer_focal_summary.csv',index=False)
    source=held[held.target_set=='Chen_held_source_library'].copy();source.to_csv(out/'source_holdout_controls.csv',index=False)
    primary=transfer[(transfer.prototype_weighting=='equal_cell_mean')&(transfer.gene_weighting=='inverse_sqrt')&(transfer.feature_mode=='all_genes')]
    controls=source.groupby('candidate_model').agg(groups=('profile','size'),matches=('matches_author_broad','sum')).reset_index()
    summary={'source_library_control':controls.to_dict('records'),
        'source_BS_groups_falsely_meiotic_with_added_LZ':int(((source.candidate_model=='Green13_ChenLZ')&(source.author_type=='BS')&(source.largest_class=='Meiotic')).sum()),
        'source_L_groups_matching_after_added_LZ':int(((source.candidate_model=='Green13_ChenLZ')&(source.author_type=='L')&source.matches_author_broad).sum()),
        'external_Leptotene_primary_largest_classes':primary[primary.author_type=='Leptotene'].largest_class.value_counts().to_dict(),
        'decision':'Do not promote any added L/Z prototype to a calibrated production reference: limited external improvement and same-source spermatogonial false positives.',
        'not_claimed':'Causal technology effects, independent donor confidence intervals, reliable in-vivo cell fractions, or resolution of older Cnot4 LZ validation.'}
    (out/'decision_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    selected=focal[(focal.target_transform=='raw')&(focal.target_weighting!='equal_cell_mean')]
    profiles=['EmptyDrops:full:Leptotene:RNA_sum','GSM5252962','GSM5252963']
    settings=[('legacy','inverse_sqrt'),('expanded','inverse_sqrt'),('legacy','uniform'),('expanded','uniform')]
    for i,p in enumerate(profiles):
        rows=selected[selected.profile==p].set_index(['gene_universe','gene_weighting'])
        axes[0].scatter(np.arange(len(settings))+(i-1)*.12,[rows.loc[key,'Meiotic'] for key in settings],label=['P15 L','Bulk L1','Bulk L2'][i])
    axes[0].set(xticks=range(4),xticklabels=['Legacy / weighted','Expanded / weighted','Legacy / uniform','Expanded / uniform'],
        ylabel='Fitted meiotic RNA weight (not cell fraction)',ylim=(0,1),title='Gene universe and weighting: external leptotene')
    axes[0].tick_params(axis='x',rotation=30);axes[0].legend()
    s=source.pivot(index=['author_type','excluded_library'],columns='candidate_model',values='Meiotic').sort_index()
    for j,model in enumerate(['Green13','Green13_ChenLZ']):axes[1].scatter(np.arange(len(s))+(j-.5)*.16,s[model],label=model)
    axes[1].set(xticks=range(len(s)),xticklabels=[f'{a}/{b}' for a,b in s.index],ylim=(0,1),
        ylabel='Fitted meiotic RNA weight (not cell fraction)',title='Held source libraries: L rescue also misassigns BS')
    axes[1].tick_params(axis='x',rotation=45);axes[1].legend()
    fig.savefig(out/'reference_limit_diagnostics.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    models=['Green13','Green13_ChenL','Green13_ChenZ','Green13_ChenLZ','Green13_ChenBS']
    for j,p in enumerate(profiles):
        rows=primary[primary.profile==p].set_index('candidate_model').reindex(models)
        for ax,col in zip(axes,['Meiotic','holdout_log2_RMSE']):ax.scatter(np.arange(len(models))+(j-1)*.13,rows[col],label=['P15 L','Bulk L1','Bulk L2'][j])
    for ax in axes:ax.set(xticks=range(len(models)),xticklabels=models);ax.tick_params(axis='x',rotation=35);ax.legend()
    axes[0].set(ylabel='Fitted meiotic RNA weight (not cell fraction)',ylim=(0,1));axes[1].set(ylabel='Held-out full-response gene log2 RMSE')
    fig.suptitle('Exploratory external transfer: equal-cell prototypes, weighted all-gene fits')
    fig.savefig(out/'candidate_transfer_diagnostics.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
