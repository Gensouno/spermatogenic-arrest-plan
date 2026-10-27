#!/usr/bin/env python3
"""Exploratory cross-study transfer of author-defined L/Z UMI prototypes."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from diagnose_reference_limits import ROOT,full_references,load_targets,normalize,extended_folds
from diagnose_reference_transfer import train_fit

MODELS={'Green13':[], 'Green13_ChenL':['L'], 'Green13_ChenZ':['Z'],
    'Green13_ChenLZ':['L','Z'], 'Green13_ChenBS':['BS']}


def component_class(column):
    if column=='Somatic':return 'Somatic'
    if column in ['GC1','Chen_BS']:return 'Spermatogonia'
    if column in [f'GC{i}' for i in range(2,9)]+['Chen_L','Chen_Z']:return 'Meiotic'
    if column in [f'GC{i}' for i in range(9,13)]:return 'Postmeiotic'
    raise ValueError('Unknown candidate reference component')


def candidate_weights(coef,columns):
    if not np.isfinite(coef).all() or (coef<0).any() or coef.sum()<=0:raise ValueError('Invalid candidate coefficients')
    p=pd.Series(coef/coef.sum(),index=columns)
    broad=pd.Series({c:sum(p[k] for k in columns if component_class(k)==c)
        for c in ['Somatic','Spermatogonia','Meiotic','Postmeiotic']})
    if not np.isclose(broad.sum(),1):raise ValueError('Broad component map lost probability mass')
    return p,broad


def specificity_panel(chen,green,per_component=200,min_cpm=10):
    """Reference-source-only selection, deterministic symbol tie-breaking."""
    chosen=set();rows=[]
    for kind,source,columns in [('Chen_stage',chen,list(chen.columns)),('Green_somatic',green,['Somatic'])]:
        abundance=normalize(source);denominator=abundance.sum(axis=1)+1e-12
        for col in columns:
            eligible=abundance[col]>=min_cpm
            table=pd.DataFrame({'gene':abundance.index[eligible],
                'specificity':(abundance[col]/denominator).loc[eligible].to_numpy(),
                'source_CPM':abundance.loc[eligible,col].to_numpy()}).sort_values(['specificity','gene'],ascending=[False,True]).head(per_component)
            for r in table.to_dict('records'):rows.append({'selection_source':kind,'component':col,**r})
            chosen.update(table.gene)
    return chosen,pd.DataFrame(rows)


def main():
    out=ROOT/'results/reference_limits';refs,old,audit=full_references();vectors,evidence,_=load_targets()
    evidence=evidence[evidence.target_set!='GSE107644'].copy()
    vectors=[v for v in vectors if v.name in set(evidence.profile)]
    if len(evidence)!=40:raise ValueError('Expected all 40 prior external profiles, with no Chen training-source profiles')
    chen=pd.read_csv(out/'chen_stage_profiles.csv.gz',index_col=0);chen_evidence=pd.read_csv(out/'chen_profile_evidence.csv')
    qc=json.loads((out/'chen_preparation_audit.json').read_text())
    if not qc['all_selected_cell_UMI_and_detected_gene_QC_match']:raise ValueError('Chen source QC not reconciled')
    genes=refs['centroid_13'].index.intersection(chen.index)
    for v in vectors:genes=genes.intersection(v.index)
    genes=genes.sort_values();foldmap=extended_folds(old['centroid_4'].index,refs['centroid_4'].index)
    folds=np.array([foldmap[g] for g in genes]);target=normalize(pd.concat([v.loc[genes] for v in vectors],axis=1))
    green=refs['centroid_13'].loc[genes];coarse=normalize(refs['centroid_4'].loc[genes])
    marker_source=chen[[f'Chen:{s}:equal_cell_mean' for s in ['A1','In','BS','BG2','G1','ePL','mPL','lPL','L','Z','eP','mP','lP','D','MI','MII','RS2','RS4','RS6','RS8']]]
    panel,panel_evidence=specificity_panel(marker_source,refs['centroid_13'])
    panel_evidence.to_csv(out/'training_marker_panel.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    panel_mask=genes.isin(panel)
    if panel_mask.sum()<500:raise ValueError('Unexpectedly small source-only feature panel')
    full_rows=[];fold_rows=[];metrics=[];condition=[];design_rows=[];prototype_rows=[]
    for view in ['equal_cell_mean','RNA_sum']:
        for stage in ['L','Z','BS']:
            name=f'Chen:{stage}:{view}';r=chen_evidence[chen_evidence.profile==name]
            if len(r)!=1 or r.iloc[0].cells<20:raise ValueError('Insufficient or ambiguous source prototype')
            prototype_rows.append({'prototype':f'Chen_{stage}','prototype_weighting':view,
                'source_cells':int(r.iloc[0].cells),'source_libraries':int(r.iloc[0].source_libraries),
                'expected_class':component_class(f'Chen_{stage}'),'source_profile':name})
        for model,stages in MODELS.items():
            reference=green.copy()
            for stage in stages:reference['Chen_'+stage]=chen.loc[genes,f'Chen:{stage}:{view}']
            reference=normalize(reference)
            reference.to_csv(out/f'candidate_{model}_{view}.csv.gz',compression={'method':'gzip','mtime':0})
            for weighting,feature_mode in [(w,f) for w in ['inverse_sqrt','uniform'] for f in ['all_genes','training_stage_markers']]:
                scale=1/np.sqrt(coarse.mean(axis=1).to_numpy()+1) if weighting=='inverse_sqrt' else np.ones(len(genes))
                fit_mask=np.ones(len(genes),dtype=bool) if feature_mode=='all_genes' else panel_mask
                s=reference.to_numpy();key={'candidate_model':model,'prototype_weighting':view,'gene_weighting':weighting,
                    'feature_mode':feature_mode,'genes':len(genes),'fit_genes':int(fit_mask.sum())}
                condition.append({**key,'components':reference.shape[1],'weighted_condition_number':float(np.linalg.cond(s[fit_mask]*scale[fit_mask,None]))})
                for row in evidence.to_dict('records'):
                    y=target[row['profile']].to_numpy();pred=np.empty_like(y);scored=pd.notna(row['expected_broad_class'])
                    for f in [-1,0,1,2,3,4]:
                        coef=train_fit(s,y,scale,fit_mask if f==-1 else (folds!=f)&fit_mask)
                        if f!=-1:pred[folds==f]=s[folds==f]@coef
                        p,broad=candidate_weights(coef,reference.columns)
                        result={**key,**row,'fold':f,'largest_class':broad.idxmax(),
                            'matches_author_broad':broad.idxmax()==row['expected_broad_class'] if scored else None,
                            'expected_class_weight':broad[row['expected_broad_class']] if scored else None,
                            **broad.to_dict(),**{f'component_{k}':float(v) for k,v in p.items()}}
                        fold_rows.append(result)
                        if f==-1:full_rows.append(result)
                    metrics.append({**key,**row,'holdout_log2_RMSE':float(np.sqrt(np.mean((np.log2(y+.5)-np.log2(pred+.5))**2))),
                        'holdout_weighted_RMSE':float(np.sqrt(np.mean(((y-pred)*scale)**2)))})
            print('Transferred',view,model,'to',len(evidence),'external profiles;',len(genes),'genes',flush=True)
    full=pd.DataFrame(full_rows);metrics=pd.DataFrame(metrics)
    full.to_csv(out/'transfer_full_weights.csv',index=False);pd.DataFrame(fold_rows).to_csv(out/'transfer_fold_weights.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    metrics.to_csv(out/'transfer_holdout_metrics.csv',index=False)
    keys=['candidate_model','prototype_weighting','gene_weighting','feature_mode','profile']
    comparison=full.merge(metrics[keys+['holdout_log2_RMSE','holdout_weighted_RMSE']],on=keys,validate='one_to_one')
    comparison.to_csv(out/'transfer_comparison.csv',index=False);pd.DataFrame(condition).to_csv(out/'transfer_conditions.csv',index=False)
    pd.DataFrame(prototype_rows).to_csv(out/'prototype_evidence.csv',index=False)
    summary=full.groupby(['candidate_model','prototype_weighting','gene_weighting','feature_mode','target_set','target_weighting']).agg(
        profiles=('profile','size'),assessable=('expected_class_weight','count'),matches=('matches_author_broad','sum')).reset_index()
    summary.to_csv(out/'transfer_label_summary.csv',index=False)
    pd.DataFrame({'gene':genes,'fold':folds}).to_csv(out/'transfer_gene_folds.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    design={'scope':'Exploratory cross-study addition of stage-defined L/Z UMI prototypes; no target-driven marker selection.',
        'target_profiles':len(evidence),'genes':len(genes),'models':MODELS,'prototype_weightings':['equal_cell_mean','RNA_sum'],
        'training_source':'Only Chen Table S2 + WT QC stage-defined L, Z and BS profiles; never P15 or external sorted labels/counts.',
        'evaluation_sources':'All 30 full P15 treatment views plus 10 prior sorted RNA profiles; no Chen source profile included in candidate evaluation.',
        'control':'Chen BS provides a same-source spermatogonial extra-column control. It is not a random negative biological population.',
        'feature_panel':{'source_only_selected_genes':len(panel),'panel_genes_in_response_universe':int(panel_mask.sum()),
            'definition':'Top 200 source-specific genes per 20 Chen stages plus 200 Green somatic-specific genes; source CPM>=10, deterministic symbol ties. Same union for all candidates; target labels/expression never rank features.',
            'evaluation':'Panel fits predict all full-response held-out genes using unchanged common-gene normalization. Panel results are a secondary follow-up, not a target-tuned final model.'},
        'gene_selection':'Deterministic shared genes across full Green references, all source prototypes and all external targets. No markers selected by target fit or class agreement.',
        'metric':'Same target values, folds and scale across candidate models within each prototype/weighting view; no held-out-specific normalization.',
        'normalization':'All common genes normalized before gene folds, so held-out totals influence column scaling; not completely prospective prediction.',
        'selection_limitation':'Motivated by previously observed external mismatch; these are cross-study exploratory transfer checks, not an untouched confirmatory set.',
        'not_claimed':['Accurate in-vivo cell-number proportions','Independent animal replication from cells or weighting views','A selected final optimal model','Causal interpretation of residuals'],
        'fixed_prepared_inputs_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [out/'chen_stage_profiles.csv.gz',out/'chen_profile_evidence.csv',out/'chen_preparation_audit.json',ROOT/'results/emptydrops_p15/type_profiles.csv.gz']},
        'generated_candidate_reference_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(out.glob('candidate_Green13*.csv.gz'))}}
    (out/'transfer_design.json').write_text(json.dumps(design,indent=2)+'\n')


if __name__=='__main__':main()
