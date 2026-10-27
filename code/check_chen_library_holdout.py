#!/usr/bin/env python3
"""Source-library holdout control and external prototype sensitivity, not animals."""
import json

import numpy as np
import pandas as pd

from diagnose_reference_limits import ROOT,full_references,load_targets,normalize,extended_folds
from check_chen_transfer import candidate_weights
from diagnose_reference_transfer import train_fit


def leave_library_profile(x,evidence,stage,excluded):
    rows=evidence[(evidence.author_type==stage)&(evidence.library!=excluded)&(evidence.weighting=='equal_cell_mean')]
    if rows.empty:raise ValueError('No remaining source library for stage')
    values=x[rows.profile].mul(rows.set_index('profile').cells,axis=1).sum(axis=1)/rows.cells.sum()
    return values,int(rows.cells.sum()),rows.library.tolist()


def main():
    out=ROOT/'results/reference_limits';refs,old,_=full_references();vectors,evidence,_=load_targets()
    evidence=evidence[evidence.target_set!='GSE107644'].copy();vectors=[v for v in vectors if v.name in set(evidence.profile)]
    batches=pd.read_csv(out/'chen_batch_prototypes.csv.gz',index_col=0);meta=pd.read_csv(out/'chen_batch_evidence.csv')
    source=meta[meta.weighting=='equal_cell_mean'];libraries=sorted(source.library.unique())
    genes=refs['centroid_13'].index.intersection(batches.index)
    for v in vectors:genes=genes.intersection(v.index)
    genes=genes.sort_values();foldmap=extended_folds(old['centroid_4'].index,refs['centroid_4'].index)
    folds=np.array([foldmap[g] for g in genes]);green=refs['centroid_13'].loc[genes]
    coarse=normalize(refs['centroid_4'].loc[genes]);scale=1/np.sqrt(coarse.mean(axis=1).to_numpy()+1)
    rows=[];conditions=[];train_audit=[];fold_rows=[]
    for excluded in libraries:
        prototypes={}
        for stage in ['L','Z']:
            v,n,remaining=leave_library_profile(batches,meta,stage,excluded);prototypes[stage]=v
            train_audit.append({'excluded_library':excluded,'stage':stage,'retained_cells':n,
                'retained_libraries':','.join(remaining),'below_20_cell_sensitivity':n<20})
        target_rows=evidence.to_dict('records');target_vectors=list(vectors)
        for r in source[source.library==excluded].itertuples():
            target_vectors.append(batches[r.profile]);target_rows.append({'profile':r.profile,'target_set':'Chen_held_source_library',
                'author_type':r.author_type,'target_weighting':'equal_cell_mean','expected_broad_class':'Spermatogonia' if r.author_type=='BS' else 'Meiotic'})
        y=normalize(pd.concat([v.loc[genes] for v in target_vectors],axis=1))
        for model in ['Green13','Green13_ChenLZ']:
            reference=green.copy()
            if model!='Green13':
                for stage in ['L','Z']:reference['Chen_'+stage]=prototypes[stage].loc[genes]
            reference=normalize(reference);s=reference.to_numpy()
            conditions.append({'excluded_library':excluded,'candidate_model':model,'genes':len(genes),
                'components':len(reference.columns),'weighted_condition_number':float(np.linalg.cond(s*scale[:,None]))})
            for r in target_rows:
                target=y[r['profile']].to_numpy();pred=np.empty_like(target);scored=pd.notna(r['expected_broad_class'])
                for f in [-1,0,1,2,3,4]:
                    coef=train_fit(s,target,scale,np.ones(len(genes),dtype=bool) if f==-1 else folds!=f)
                    if f!=-1:pred[folds==f]=s[folds==f]@coef
                    p,broad=candidate_weights(coef,reference.columns)
                    result={'excluded_library':excluded,'candidate_model':model,**r,'fold':f,'largest_class':broad.idxmax(),
                        'matches_author_broad':broad.idxmax()==r['expected_broad_class'] if scored else None,
                        'expected_class_weight':broad[r['expected_broad_class']] if scored else None,**broad.to_dict(),
                        **{f'component_{k}':float(v) for k,v in p.items()}}
                    fold_rows.append(result)
                    if f==-1:row=result.copy()
                rows.append({**row,'holdout_log2_RMSE':float(np.sqrt(np.mean((np.log2(target+.5)-np.log2(pred+.5))**2))),
                    'holdout_weighted_RMSE':float(np.sqrt(np.mean(((target-pred)*scale)**2)))})
        print('Checked leave-source-library-out',excluded,len(target_rows),'targets',flush=True)
    table=pd.DataFrame(rows);table.to_csv(out/'source_library_holdout.csv',index=False)
    pd.DataFrame(fold_rows).to_csv(out/'source_library_holdout_folds.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    pd.DataFrame(conditions).to_csv(out/'source_library_conditions.csv',index=False)
    pd.DataFrame(train_audit).to_csv(out/'source_library_training_inventory.csv',index=False)
    external=table[table.target_set!='Chen_held_source_library']
    external.groupby(['candidate_model','profile','target_set','author_type','target_weighting']).agg(
        library_variants=('excluded_library','size'),Meiotic_min=('Meiotic','min'),Meiotic_max=('Meiotic','max'),
        RMSE_min=('holdout_log2_RMSE','min'),RMSE_max=('holdout_log2_RMSE','max'),
        largest_class_variants=('largest_class','nunique')).reset_index().to_csv(out/'source_library_external_ranges.csv',index=False)
    design={'scope':'Leave one source processing library out of both Chen L/Z prototypes; unmodified equal-cell means, all fixed shared genes and inverse-sqrt Green4 weights.',
        'source_libraries':libraries,'held_source_stage_groups':len(source),'external_profiles':len(evidence),'genes':len(genes),
        'source_positive_control':'Held-library stage pseudobulks contain no cells used in that variant L/Z prototype. Broad labels from author S2. This controls source-library generalization only, not independent animals or cross-study calibration.',
        'low_N':'Retained source groups below 20 cells explicitly marked as sensitivity, not a selected production reference.',
        'no_panel_leakage':'No source-selected marker panel in this stage. Full fixed gene universe and independent Green weights only.',
        'interpretation':'External variants share the same source and target cells; ranges are processing sensitivity, not confidence intervals. Same-source success does not validate external reference transfer or cell-number fractions.'}
    (out/'source_library_design.json').write_text(json.dumps(design,indent=2)+'\n')


if __name__=='__main__':main()
