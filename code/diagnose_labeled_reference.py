#!/usr/bin/env python3
"""Fixed-reference prediction against independently supplied author cell-type labels."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from diagnose_reference_transfer import MODELS,broad_weights,gene_folds,train_fit

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/labeled_reference';out.mkdir(parents=True,exist_ok=True)
    x=pd.read_csv(out/'type_profiles.csv.gz',index_col=0)
    evidence=pd.read_csv(out/'profile_evidence.csv')
    if x.columns.tolist()!=evidence.profile.tolist() or x.index.duplicated().any():raise ValueError('Profile inventory differs')
    refs={m:pd.read_csv(ROOT/f'results/stage_reference/pooled_{m}_reference.csv.gz',index_col=0) for m in MODELS}
    foldmap=gene_folds(sorted(refs['centroid_4'].index))
    genes=refs['centroid_4'].index.intersection(x.index).sort_values()
    y=x.loc[genes].div(x.loc[genes].sum(axis=0),axis=1)*1e6
    refs={m:r.loc[genes].div(r.loc[genes].sum(axis=0),axis=1)*1e6 for m,r in refs.items()}
    scale=1/np.sqrt(refs['centroid_4'].mean(axis=1).to_numpy()+1)
    folds=np.array([foldmap[g] for g in genes]);metrics=[];weights=[];conditions=[]
    for model,ref in refs.items():
        s=ref.to_numpy()
        conditions.append({'model':model,'genes':len(genes),'weighted_condition_number':float(np.linalg.cond(s*scale[:,None]))})
        for row in evidence.to_dict('records'):
            target=y[row['profile']].to_numpy();prediction=np.empty_like(target)
            for fold in [-1,0,1,2,3,4]:
                coef=train_fit(s,target,scale,np.ones(len(genes),dtype=bool) if fold==-1 else folds!=fold)
                if fold!=-1:prediction[folds==fold]=s[folds==fold]@coef
                p=broad_weights(coef,ref.columns,model)
                weights.append({**row,'model':model,'fold':fold,'largest_class':p.idxmax(),
                    'matches_author_broad':p.idxmax()==row['expected_broad_class'],
                    'expected_class_weight':p[row['expected_broad_class']],**p.to_dict()})
            metrics.append({**row,'model':model,'genes':len(genes),
                'shared_gene_expression_mass':float(x.loc[genes,row['profile']].sum()/x[row['profile']].sum()),
                'holdout_log2_RMSE':float(np.sqrt(np.mean((np.log2(target+.5)-np.log2(prediction+.5))**2))),
                'holdout_weighted_RMSE':float(np.sqrt(np.mean(((target-prediction)*scale)**2)))})
        print('Diagnosed',model,len(evidence),'author-labelled profiles;',len(genes),'genes',flush=True)
    metrics=pd.DataFrame(metrics);weights=pd.DataFrame(weights)
    metrics.to_csv(out/'holdout_metrics.csv',index=False)
    weights.to_csv(out/'fold_weights.csv',index=False)
    full=weights[weights.fold==-1].copy();full.to_csv(out/'full_fit_weights.csv',index=False)
    summary=full.groupby(['Sample','Library','weighting','model']).agg(
        type_profiles=('profile','size'),matches=('matches_author_broad','sum'),
        median_author_class_weight=('expected_class_weight','median')).reset_index()
    summary=summary.merge(metrics.groupby(['Sample','Library','weighting','model']).holdout_log2_RMSE.median().reset_index(),
        on=['Sample','Library','weighting','model'],validate='one_to_one')
    summary.to_csv(out/'library_summary.csv',index=False)
    full.groupby(['Sample','AnnotatedClusters','expected_broad_class','weighting','model']).agg(
        matches=('matches_author_broad','sum'),profiles=('profile','size'),
        median_author_class_weight=('expected_class_weight','median')).to_csv(out/'label_summary.csv')
    weights[weights.fold!=-1].groupby(['Sample','weighting','model','fold']).agg(
        type_profiles=('profile','size'),matches=('matches_author_broad','sum')).to_csv(out/'fold_summary.csv')
    pd.DataFrame(conditions).to_csv(out/'design_conditions.csv',index=False)
    inputs=[f'results/stage_reference/pooled_{m}_reference.csv.gz' for m in MODELS]
    notes={'scope':'Three preselected complete normal libraries: adult B6 do17815, P5 do26386, P15 do18195; author-labelled library/type groups with >=20 cells; two expression weighting views, not independent repeats.',
        'label_validation':'Author transcriptomic annotations are external to this project fitting, not independently measured histological purity or cell-number ground truth.',
        'reference_training':'No new reference trained and no markers retuned. Frozen prior 3 references and gene universe reused.',
        'gene_folds_seed':20261024,'genes':len(genes),
        'normalization':'Shared-gene target/reference columns normalized once before gene splitting; held-out genes therefore contribute to column totals. Not a completely prospective missing-gene experiment.',
        'prediction':'Four folds train, fifth predicts, no held-out-specific rescaling. Same weights and folds as prior transfer diagnosis.',
        'not_claimed':['Complete atlas ingestion','Independent animal replication','Accurate in-vivo cell-number fractions','Age or genotype causal effects','New calibrated age-matched reference'],
        'fixed_inputs_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in inputs}}
    (out/'design.json').write_text(json.dumps(notes,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    subset=full[full.weighting=='RNA_sum']
    fig,axes=plt.subplots(1,3,figsize=(13,5),layout='constrained')
    for ax,sample in zip(axes,['P5','P15','B6']):
        table=subset[subset.Sample==sample].pivot(index='AnnotatedClusters',columns='model',values='expected_class_weight').reindex(columns=MODELS)
        image=ax.imshow(table.to_numpy(),aspect='auto',vmin=0,vmax=1,cmap='viridis')
        ax.set(xticks=range(3),xticklabels=['Original 4','Centroid 4','Centroid 13'],
            yticks=range(len(table)),yticklabels=table.index,title=sample+' author-labelled groups')
        ax.tick_params(axis='x',rotation=45)
    fig.colorbar(image,ax=axes,label='Fitted author broad-class weight (not cell fraction)',shrink=.8)
    fig.suptitle('External author-label compatibility: three selected libraries, RNA-sum profiles')
    fig.savefig(out/'author_label_compatibility.png',dpi=180);plt.close(fig)
    snap=ROOT/'labeled_reference_snapshot';snap.mkdir(exist_ok=True)
    names=['all_group_inventory.csv','profile_evidence.csv','library_evidence.csv','preparation_audit.json',
        'holdout_metrics.csv','full_fit_weights.csv','library_summary.csv','label_summary.csv','fold_summary.csv',
        'design_conditions.csv','design.json','author_label_compatibility.png']
    names.append('emptydrops_label_audit.json')
    names+=['zygotene_overlap_audit.json','zygotene_overlap_metrics.csv','zygotene_overlap_weights.csv']
    for name in names:shutil.copyfile(out/name,snap/name)
    sources=['scripts/fetch_labeled_reference.py','scripts/fetch_labeled_ranges.py','scripts/run_labeled_reference.py',
        'scripts/audit_emptydrops_labels.py',
        'scripts/check_labeled_reproduction.py','scripts/check_labeled_snapshot.py','tests/test_labeled_reference.py',
        'code/prepare_labeled_reference.py','code/diagnose_labeled_reference.py','code/diagnose_reference_transfer.py',
        'code/check_author_zygotene.py',
        'code/check_sorted_expansion.py','scripts/discover_expansion.py','scripts/check_project.py',
        'data/labeled_reference_contract.json']
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest={'scope':notes['scope'],'source_sha256':{n:digest(ROOT/n) for n in sources},
        'input_specification_sha256':{n:digest(ROOT/n) for n in ['data/labeled_reference_assets.json','data/emptydrops_label_assets.json']},
        'files':[{'path':str(p.relative_to(snap)),'bytes':p.stat().st_size,'sha256':digest(p)}
            for p in sorted(snap.rglob('*')) if p.is_file() and p.name!='manifest.json']}
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(summary.to_string(index=False))


if __name__=='__main__':main()
