#!/usr/bin/env python3
"""Separate recheck of all previously audited sorted profiles, including Cnot4 LZ."""
import json

import numpy as np
import pandas as pd

from diagnose_reference_limits import ROOT,full_references,normalize,extended_folds
from check_sorted_expansion import tar_profiles
from check_chen_transfer import candidate_weights,MODELS
from diagnose_reference_transfer import train_fit
import sys
sys.path.insert(0,str(ROOT/'scripts'))
from freeze_reference_limits_inputs import digest


def main():
    out=ROOT/'results/reference_limits';old_design=ROOT/'reference_audit_snapshot/design.json'
    prior=json.loads(old_design.read_text())['fixed_inputs_sha256']
    names=[k for k in prior if any(t in k for t in ['GSE158471','GSE232415','H3mm15_sorted_symbol_counts','source_sample_evidence'])]
    if len(names)!=6:raise ValueError('Expected six old sorted source/prepared inputs')
    assets=[]
    for name in names:
        p=ROOT/name
        if digest(p)!=prior[name]:raise ValueError('Prior sorted input differs')
        assets.append({'path':name,'bytes':p.stat().st_size,'sha256':prior[name]})
    assets.append({'path':str(old_design.relative_to(ROOT)),'bytes':old_design.stat().st_size,'sha256':digest(old_design)})
    spec=ROOT/'data/reference_limits_sorted_inputs.json';value={'scope':'Previously cached sorted-profile inputs checked against the pre-existing reference audit; no new network retrieval.','assets':assets}
    if spec.exists() and json.loads(spec.read_text())!=value:raise ValueError('Frozen sorted source input specification changed')
    if not spec.exists():spec.write_text(json.dumps(value,indent=2)+'\n')
    vectors=[];rows=[]
    for gse,selected,unit in [('GSE158471',{'GSM4801338':'LZ','GSM4801339':'LZ','GSM4801340':'PD','GSM4801341':'PD'},'FPKM'),
        ('GSE232415',{'GSM8618769':'PS','GSM8618770':'PS','GSM8618771':'RS','GSM8618772':'RS'},'counts')]:
        x,meta=tar_profiles(gse,selected,unit)
        for r in meta:
            vectors.append(x[r['gsm']].rename(r['gsm']));rows.append({'profile':r['gsm'],'gse':gse,'stage':r['stage'],'unit':unit,
                'expected_broad_class':'Postmeiotic' if r['stage']=='RS' else 'Meiotic'})
    x=pd.read_csv(ROOT/'data/interim/expansion/H3mm15_sorted_symbol_counts.tsv.gz',sep='\t',index_col=0)
    e=pd.read_csv(ROOT/'results/data_expansion/sorted_validation/source_sample_evidence.csv');e=e[e.gse=='GSE290629']
    if len(e)!=12:raise ValueError('Expected 12 frozen H3mm15 source profiles')
    for r in e.itertuples():
        vectors.append(x[r.original_column].rename(r.gsm));rows.append({'profile':r.gsm,'gse':'GSE290629','stage':r.stage,'unit':'counts',
            'expected_broad_class':None if r.stage=='preL' else ('Spermatogonia' if r.stage=='SG' else 'Postmeiotic' if r.stage=='RS' else 'Meiotic')})
    evidence=pd.DataFrame(rows)
    references={m:pd.read_csv(out/f'candidate_{m}_equal_cell_mean.csv.gz',index_col=0) for m in MODELS}
    genes=references['Green13'].index
    for v in vectors:genes=genes.intersection(v.index)
    genes=genes.sort_values();y=normalize(pd.concat([v.loc[genes] for v in vectors],axis=1))
    full,old,_=full_references();coarse=normalize(full['centroid_4'].loc[genes])
    foldmap=extended_folds(old['centroid_4'].index,full['centroid_4'].index);folds=np.array([foldmap[g] for g in genes])
    panel=set(pd.read_csv(out/'training_marker_panel.csv.gz').gene);panel_mask=genes.isin(panel)
    records=[];fold_rows=[];conditions=[]
    for weighting in ['inverse_sqrt','uniform']:
        scale=1/np.sqrt(coarse.mean(axis=1).to_numpy()+1) if weighting=='inverse_sqrt' else np.ones(len(genes))
        for feature in ['all_genes','training_stage_markers']:
            mask=np.ones(len(genes),dtype=bool) if feature=='all_genes' else panel_mask
            for model,reference in references.items():
                ref=normalize(reference.loc[genes]);s=ref.to_numpy()
                key={'candidate_model':model,'gene_weighting':weighting,'feature_mode':feature,'genes':len(genes),'fit_genes':int(mask.sum())}
                conditions.append({**key,'weighted_condition_number':float(np.linalg.cond(s[mask]*scale[mask,None]))})
                for row in rows:
                    t=y[row['profile']].to_numpy();prediction=np.empty_like(t);scored=row['expected_broad_class'] is not None
                    for f in [-1,0,1,2,3,4]:
                        coef=train_fit(s,t,scale,mask if f==-1 else mask&(folds!=f))
                        if f!=-1:prediction[folds==f]=s[folds==f]@coef
                        p,b=candidate_weights(coef,ref.columns)
                        result={**key,**row,'fold':f,'largest_class':b.idxmax(),
                            'matches_author_broad':b.idxmax()==row['expected_broad_class'] if scored else None,
                            'expected_class_weight':b[row['expected_broad_class']] if scored else None,**b.to_dict()}
                        fold_rows.append(result)
                        if f==-1:record=result.copy()
                    vector=next(v for v in vectors if v.name==row['profile'])
                    records.append({**record,'shared_gene_expression_mass':float(vector.loc[genes].sum()/vector.sum()),
                        'holdout_log2_RMSE':float(np.sqrt(np.mean((np.log2(t+.5)-np.log2(prediction+.5))**2)))})
    pd.DataFrame(records).to_csv(out/'legacy_sorted_comparison.csv',index=False)
    pd.DataFrame(fold_rows).to_csv(out/'legacy_sorted_fold_weights.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    pd.DataFrame(conditions).to_csv(out/'legacy_sorted_conditions.csv',index=False);evidence.to_csv(out/'legacy_sorted_samples.csv',index=False)
    pd.DataFrame({'gene':genes,'fold':folds}).to_csv(out/'legacy_sorted_gene_design.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    design={'scope':'Separate exploratory recheck of all 20 prior sorted WT profiles, two preL unscored; no untouched-validation claim.',
        'profiles':len(rows),'assessable':18,'genes':len(genes),'source_only_panel_genes':int(panel_mask.sum()),
        'candidate_prototypes':'Fixed equal-cell means only; same five candidates as the primary transfer branch.',
        'normalization':'Common-gene response normalized before gene folds; within-branch candidates have identical response, genes and scales. Do not directly compare cross-branch RMSE with different gene universes.',
        'source_status':'Three study sources, mixed LZ not pure L, FPKM vs counts and purification/age confounders retained.',
        'input_hashes_verified_against':'reference_audit_snapshot/design.json',
        'fixed_candidate_reference_sha256':{str(p.relative_to(ROOT)):digest(p) for p in sorted(out.glob('candidate_Green13*_equal_cell_mean.csv.gz'))}}
    (out/'legacy_sorted_design.json').write_text(json.dumps(design,indent=2)+'\n')
    print('Rechecked all',len(rows),'old sorted profiles;',len(genes),'common genes',flush=True)


if __name__=='__main__':main()
