#!/usr/bin/env python3
"""Fixed reference comparison across complete, new and matched P15 cells."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from diagnose_reference_transfer import MODELS,broad_weights,gene_folds,train_fit

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/emptydrops_p15'
    x=pd.read_csv(out/'type_profiles.csv.gz',index_col=0);evidence=pd.read_csv(out/'profile_evidence.csv')
    if x.columns.tolist()!=evidence.profile.tolist() or x.index.duplicated().any():raise ValueError('Profile identity differs')
    paths=[f'results/stage_reference/pooled_{m}_reference.csv.gz' for m in MODELS]
    refs={m:pd.read_csv(ROOT/p,index_col=0) for m,p in zip(MODELS,paths)}
    foldmap=gene_folds(sorted(refs['centroid_4'].index))
    genes=refs['centroid_4'].index.intersection(x.index).sort_values()
    if len(genes)<10000:raise ValueError('Unexpectedly low fixed-reference gene coverage')
    y=x.loc[genes].div(x.loc[genes].sum(),axis=1)*1e6
    refs={m:r.loc[genes].div(r.loc[genes].sum(),axis=1)*1e6 for m,r in refs.items()}
    scale=1/np.sqrt(refs['centroid_4'].mean(axis=1).to_numpy()+1)
    folds=np.array([foldmap[g] for g in genes]);weights=[];metrics=[];conditions=[]
    for model,ref in refs.items():
        s=ref.to_numpy();conditions.append({'model':model,'genes':len(genes),'weighted_condition_number':float(np.linalg.cond(s*scale[:,None]))})
        for row in evidence.to_dict('records'):
            target=y[row['profile']].to_numpy();prediction=np.empty_like(target)
            shared_mass=float(x.loc[genes,row['profile']].sum()/x[row['profile']].sum())
            for fold in [-1,0,1,2,3,4]:
                train=np.ones(len(genes),dtype=bool) if fold==-1 else folds!=fold
                coef=train_fit(s,target,scale,train)
                if fold!=-1:prediction[folds==fold]=s[folds==fold]@coef
                p=broad_weights(coef,ref.columns,model)
                weights.append({**row,'model':model,'fold':fold,'largest_class':p.idxmax(),
                    'matches_author_broad':p.idxmax()==row['expected_broad_class'],
                    'expected_class_weight':p[row['expected_broad_class']],**p.to_dict()})
            metrics.append({**row,'model':model,'genes':len(genes),
                'shared_gene_expression_mass':shared_mass,
                'reference_gene_expression_mass_of_unique_mouse':shared_mass*row['branch_common_expression_mass'],
                'reference_gene_expression_mass_of_all_features':shared_mass*row['branch_common_expression_mass']*row['unique_mouse_expression_mass'],
                'holdout_log2_RMSE':float(np.sqrt(np.mean((np.log2(target+.5)-np.log2(prediction+.5))**2))),
                'holdout_weighted_RMSE':float(np.sqrt(np.mean(((target-prediction)*scale)**2)))})
        print('Diagnosed',model,len(evidence),'profiles;',len(genes),'shared genes',flush=True)
    metrics=pd.DataFrame(metrics);weights=pd.DataFrame(weights)
    metrics.to_csv(out/'holdout_metrics.csv',index=False);weights.to_csv(out/'fold_weights.csv',index=False)
    full=weights[weights.fold==-1].copy();full.to_csv(out/'full_fit_weights.csv',index=False)
    keys=['branch','stratum','weighting','model']
    summary=full.groupby(keys).agg(type_profiles=('profile','size'),matches=('matches_author_broad','sum'),
        median_author_class_weight=('expected_class_weight','median')).reset_index().merge(
        metrics.groupby(keys).holdout_log2_RMSE.median().reset_index(),on=keys,validate='one_to_one')
    summary.to_csv(out/'stratum_summary.csv',index=False)
    focal=full[full.author_type.isin(['Leptotene','Zygotene','Pachytene_1','Pachytene_2'])].merge(
        metrics[['profile','model','holdout_log2_RMSE','holdout_weighted_RMSE']],on=['profile','model'],validate='one_to_one')
    focal.to_csv(out/'early_meiotic_comparison.csv',index=False)
    pairing=['stratum','author_type','weighting','model']
    paired=full[full.branch=='EmptyDrops'].merge(full[full.branch=='CellRanger'],on=pairing,suffixes=('_ED','_CR'),validate='one_to_one')
    for key in ['expected_class_weight','cells']:
        if key=='cells' and not paired.cells_ED.eq(paired.cells_CR).all():raise ValueError('Matched cells differ')
    paired['expected_weight_ED_minus_CR']=paired.expected_class_weight_ED-paired.expected_class_weight_CR
    paired=paired.merge(metrics[['profile','model','holdout_log2_RMSE']].rename(columns={'profile':'profile_ED','holdout_log2_RMSE':'RMSE_ED'}),on=['profile_ED','model'])
    paired=paired.merge(metrics[['profile','model','holdout_log2_RMSE']].rename(columns={'profile':'profile_CR','holdout_log2_RMSE':'RMSE_CR'}),on=['profile_CR','model'])
    paired['RMSE_ED_minus_CR']=paired.RMSE_ED-paired.RMSE_CR
    paired.to_csv(out/'matched_count_branch_comparison.csv',index=False)
    pd.DataFrame(conditions).to_csv(out/'design_conditions.csv',index=False)
    pd.DataFrame({'gene':genes,'fold':folds,'weight':scale}).to_csv(out/'gene_design.csv.gz',index=False)
    notes={'scope':'One complete P15 EmptyDrops library; all eligible author groups, complete/new/overlap strata, two weighting views and same-barcode CellRanger comparison.',
        'genes':len(genes),'gene_folds_seed':20261024,'labels':'Original author P15Clusters, not independently measured purity or cell-number truth.',
        'normalization':'Identical common genes for both count branches and every profile. Column totals normalized once before gene splitting; held-out genes contribute to totals.',
        'references':'Three unchanged earlier references, fixed inverse-square-root centroid-4 weights, fixed five gene folds; no marker retuning or new reference training.',
        'interpretation':'RNA coefficients are compatibility diagnostics, not in-vivo cell-number fractions. Largest-class matches are not calibrated classification accuracy.',
        'selection':'Follow-up prompted by earlier overlapping Zygotene result, not untouched independent validation.',
        'not_claimed':['Full atlas ingestion','Biological replication','Causal age/genotype effects','Reliable L/Z fractions','Resolution of Cnot4 LZ mismatch','Prospective missing-gene prediction'],
        'fixed_inputs_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}
    (out/'design.json').write_text(json.dumps(notes,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    colors={'original_4':'#777777','centroid_4':'#d58223','centroid_13':'#3279a8'}
    selected=focal[(focal.branch=='EmptyDrops')&(focal.weighting=='RNA_sum')&focal.author_type.isin(['Leptotene','Zygotene'])]
    names=[f'{t}/{s}' for t in ['Leptotene','Zygotene'] for s in ['full','overlap','EmptyDrops_only']
        if ((selected.author_type==t)&(selected.stratum==s)).any()]
    for model in MODELS:
        table=selected[selected.model==model].copy();table['name']=table.author_type+'/'+table.stratum;table=table.set_index('name').reindex(names)
        for ax,col in zip(axes,['expected_class_weight','holdout_log2_RMSE']):
            ax.plot(range(len(names)),table[col],marker='o',label=model,color=colors[model])
    for ax in axes:ax.set(xticks=range(len(names)),xticklabels=names);ax.tick_params(axis='x',rotation=40);ax.legend()
    axes[0].set(ylabel='Fitted meiotic RNA weight (not cell fraction)',ylim=(0,1))
    axes[1].set(ylabel='Held-out gene log2 RMSE')
    fig.suptitle('Complete P15 early-meiotic author labels: RNA-sum profiles, one library')
    fig.savefig(out/'early_meiotic_diagnostics.png',dpi=170);plt.close(fig)
    print(summary.to_string(index=False),flush=True)


if __name__=='__main__':main()
