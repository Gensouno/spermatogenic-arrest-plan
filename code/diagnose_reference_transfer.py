#!/usr/bin/env python3
"""Fixed-reference gene-holdout diagnostics; never proportion-accuracy claims."""
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from check_sorted_expansion import tar_profiles

ROOT=Path(__file__).resolve().parents[1]
MODELS=['original_4','centroid_4','centroid_13']
EXPECTED={'SG':'Spermatogonia','LZ':'Meiotic','PD':'Meiotic','PS':'Meiotic','Mll':'Meiotic','RS':'Postmeiotic'}


def gene_folds(genes):
    if len(genes)!=len(set(genes)) or len(genes)<5:raise ValueError('Need unique genes and at least five features')
    return dict(zip(genes,np.random.default_rng(20261024).permutation(len(genes))%5))


def train_fit(reference,target,scale,train):
    if (reference.ndim!=2 or len(target)!=len(reference) or len(scale)!=len(target)
            or train.shape!=target.shape or train.dtype!=bool or not train.any()
            or not np.isfinite(reference).all() or not np.isfinite(target).all()
            or (reference<0).any() or (target<0).any() or not np.isfinite(scale).all() or (scale<=0).any()):
        raise ValueError('Invalid fixed-reference fit')
    coef,_=nnls(reference[train]*scale[train,None],target[train]*scale[train],maxiter=2000)
    if coef.sum()<=0:raise ValueError('Empty fit')
    return coef


def broad_weights(coef,columns,model):
    p=pd.Series(coef/coef.sum(),index=columns)
    if model!='centroid_13':return p
    return pd.Series({'Somatic':p['Somatic'],'Spermatogonia':p['GC1'],
        'Meiotic':p[[f'GC{i}' for i in range(2,9)]].sum(),
        'Postmeiotic':p[[f'GC{i}' for i in range(9,13)]].sum()})


def main():
    out=ROOT/'results/reference_transfer';out.mkdir(parents=True,exist_ok=True)
    fixed={m:pd.read_csv(ROOT/f'results/stage_reference/pooled_{m}_reference.csv.gz',index_col=0) for m in MODELS}
    universe=sorted(fixed['centroid_4'].index);foldmap=gene_folds(universe)
    pd.DataFrame({'gene':universe,'fold':[foldmap[g] for g in universe]}).to_csv(out/'gene_folds.csv',index=False)
    targets=[];inputs={str((ROOT/f'results/stage_reference/pooled_{m}_reference.csv.gz').relative_to(ROOT)):
        hashlib.sha256((ROOT/f'results/stage_reference/pooled_{m}_reference.csv.gz').read_bytes()).hexdigest() for m in MODELS}
    for gse,selected,unit in [
        ('GSE158471',{'GSM4801338':'LZ','GSM4801339':'LZ','GSM4801340':'PD','GSM4801341':'PD'},'FPKM'),
        ('GSE232415',{'GSM8618769':'PS','GSM8618770':'PS','GSM8618771':'RS','GSM8618772':'RS'},'counts')]:
        x,meta=tar_profiles(gse,selected,unit)
        targets.append((gse,x,[{'profile':r['gsm'],'stage':r['stage'],'kind':'sorted_WT'} for r in meta]))
        for p in [ROOT/f'data/raw/expansion_assets/{gse}/{gse}_RAW.tar',ROOT/f'results/data_expansion/{gse}_metadata.json']:
            inputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    x=pd.read_csv(ROOT/'data/interim/expansion/H3mm15_sorted_symbol_counts.tsv.gz',sep='\t',index_col=0)
    evidence=pd.read_csv(ROOT/'results/data_expansion/sorted_validation/source_sample_evidence.csv')
    meta=evidence[evidence.gse=='GSE290629']
    x=x[meta.original_column.tolist()];x.columns=meta.gsm.tolist()
    targets.append(('GSE290629',x,[{'profile':r.gsm,'stage':r.stage,'kind':'sorted_WT'} for _,r in meta.iterrows()]))
    for p in [ROOT/'data/interim/expansion/H3mm15_sorted_symbol_counts.tsv.gz',
              ROOT/'results/data_expansion/sorted_validation/source_sample_evidence.csv']:
        inputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    age=pd.read_csv(ROOT/'results/age_reference/library_profiles.csv.gz',index_col=0)
    targets.append(('GSE121904',age,[{'profile':c,'stage':'unlabeled_capture','kind':'age_preparation_capture'} for c in age]))
    metrics=[];weights=[];conditions=[]
    for dataset,x,meta in targets:
        genes=fixed['centroid_4'].index.intersection(x.index).sort_values()
        refs={m:r.loc[genes].div(r.loc[genes].sum(),axis=1)*1e6 for m,r in fixed.items()}
        y=x.loc[genes].div(x.loc[genes].sum(),axis=1)*1e6
        scale=1/np.sqrt(refs['centroid_4'].mean(axis=1).to_numpy()+1)
        folds=np.array([foldmap[g] for g in genes])
        for model,ref in refs.items():
            s=ref.to_numpy();condition=float(np.linalg.cond(s*scale[:,None]))
            conditions.append({'dataset':dataset,'model':model,'genes':len(genes),'weighted_condition_number':condition})
            for row in meta:
                target=y[row['profile']].to_numpy();prediction=np.empty_like(target)
                full=train_fit(s,target,scale,np.ones(len(genes),dtype=bool))
                for fold in [-1,0,1,2,3,4]:
                    coef=full if fold==-1 else train_fit(s,target,scale,folds!=fold)
                    if fold!=-1:prediction[folds==fold]=s[folds==fold]@coef
                    p=broad_weights(coef,ref.columns,model)
                    expected=EXPECTED.get(row['stage']) if row['kind']=='sorted_WT' else None
                    weights.append({'dataset':dataset,**row,'model':model,'fold':fold,
                        'largest_class':p.idxmax(),'expected_broad_class':expected,
                        'expected_class_weight':p.get(expected,np.nan),**p.to_dict()})
                metrics.append({'dataset':dataset,**row,'model':model,'genes':len(genes),
                    'holdout_log2_RMSE':float(np.sqrt(np.mean((np.log2(target+.5)-np.log2(prediction+.5))**2))),
                    'holdout_weighted_RMSE':float(np.sqrt(np.mean(((target-prediction)*scale)**2))),
                    'training_weighted_RMSE':float(np.sqrt(np.mean(((target-s@full)*scale)**2)))})
        print('Finished',dataset,len(meta),'profiles;',len(genes),'shared genes',flush=True)
    metrics=pd.DataFrame(metrics);weights=pd.DataFrame(weights)
    metrics.to_csv(out/'holdout_metrics.csv',index=False);weights.to_csv(out/'fold_weights.csv',index=False)
    pd.DataFrame(conditions).to_csv(out/'design_conditions.csv',index=False)
    scored=weights[(weights.kind=='sorted_WT')&weights.expected_broad_class.notna()].copy()
    scored['matches']=scored.largest_class==scored.expected_broad_class
    scored.groupby(['dataset','model','fold']).agg(profiles=('profile','size'),matches=('matches','sum')).to_csv(out/'sorted_class_matches.csv')
    ranges=weights.groupby(['dataset','profile','model'])[['Somatic','Spermatogonia','Meiotic','Postmeiotic']].agg(lambda v:float(v.max()-v.min()))
    ranges.to_csv(out/'weight_ranges.csv')
    comparison=metrics.pivot(index=['dataset','profile','kind'],columns='model',values='holdout_log2_RMSE')
    comparison['fine_minus_broad']=comparison.centroid_13-comparison.centroid_4
    comparison.to_csv(out/'fine_vs_broad.csv')
    notes={'scope':'Fixed reference compatibility on 20 sorted WT profiles and 32 age/preparation/QC/weighting profiles; not 52 independent animals.',
        'gene_folds_seed':20261024,'normalization':'Each shared-gene target/reference column is normalized once to 1e6 before gene splitting; total normalization uses all shared genes, not a fully prospective missing-gene experiment.',
        'predictions':'Coefficients fit on four folds; fifth predicted with no held-out-gene-specific calibration or rescaling.',
        'gene_universe':'Frozen prior reference universe selected through availability in an earlier external table; no target marker selection.',
        'limitations':['Captured libraries lack independently supplied cell-type truth labels.',
            'Age, enrichment and preparation confounded; no age p-values or cell-number fractions.',
            'Holding genes out checks predictive compatibility, not proportion accuracy.',
            'The H3mm15 gene-symbol matrix is a frozen previously prepared input, not reannotated in this run.',
            'Fold weight ranges are stability diagnostics, not confidence intervals.'],
        'fixed_inputs_sha256':inputs}
    (out/'design.json').write_text(json.dumps(notes,indent=2)+'\n')
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    order=['GSE158471','GSE232415','GSE290629','GSE121904']
    for j,m in enumerate(MODELS):
        vals=[metrics[(metrics.dataset==d)&(metrics.model==m)].holdout_log2_RMSE.median() for d in order]
        ax.bar(np.arange(4)+(j-1)*.24,vals,width=.24,label=m)
    ax.set(xticks=range(4),xticklabels=['Cnot4 sorted','NKAPL sorted','H3mm15 sorted','Age/preparation captures'],
           ylabel='Median gene-holdout log2 RMSE',title='Fixed-reference predictive compatibility (not proportion accuracy)')
    ax.legend(frameon=False);fig.savefig(out/'holdout_summary.png',dpi=180);plt.close(fig)
    snap=ROOT/'reference_audit_snapshot';snap.mkdir(exist_ok=True)
    for name in ['holdout_metrics.csv','fold_weights.csv','design_conditions.csv','sorted_class_matches.csv',
                 'weight_ranges.csv','fine_vs_broad.csv','design.json','holdout_summary.png']:
        shutil.copyfile(out/name,snap/name)
    for name in ['source_audit.json','library_QC.csv','author_sample_table.csv','profile_evidence.csv']:
        shutil.copyfile(ROOT/'results/age_reference'/name,snap/name)
    sources=['scripts/fetch_age_reference.py','scripts/run_reference_audit.py','scripts/check_reference_reproduction.py',
             'scripts/check_project.py','tests/test_reference_transfer.py','code/prepare_age_reference.py','code/diagnose_reference_transfer.py',
             'code/check_sorted_expansion.py','scripts/discover_expansion.py','data/reference_audit_contract.json']
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest={'scope':notes['scope'],'source_sha256':{n:digest(ROOT/n) for n in sources},
        'input_specification_sha256':{'data/age_reference_assets.json':digest(ROOT/'data/age_reference_assets.json')},
        'files':[{'path':str(p.relative_to(snap)),'bytes':p.stat().st_size,'sha256':digest(p)}
                 for p in sorted(snap.rglob('*')) if p.is_file() and p.name!='manifest.json']}
    (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(comparison.groupby('dataset').median().to_string())


if __name__=='__main__':main()
