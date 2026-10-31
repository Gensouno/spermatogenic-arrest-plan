#!/usr/bin/env python3
"""Outcome-blind reference neighborhoods for frozen descriptive priorities."""
import argparse
import json
from pathlib import Path
import platform
import importlib.metadata
import numpy as np
import pandas as pd
from check_project import verify_manifest
from run_independent_contexts import prepare,digest,ROOT,COHORTS

REFERENCES={'Chen':{'Spermatogonia','Preleptotene','Meiotic','Postmeiotic'},
            'Green':{'Spermatogonia','Meiotic_or_transition','Postmeiotic'}}
COMPARISONS=['abundance_only','stage_abundance_dominance']


def stage_available(frame,reference):
    return (frame[reference+'_measured'].eq(True)&frame[reference+'_peak_group'].isin(REFERENCES[reference])
            &frame[reference+'_peak_group_mean_CPM'].gt(0)
            &np.isfinite(frame[reference+'_group_dominance_share']))


def neighbors(target,background,reference,caliper,strict):
    """Never use expression effects or perturbation abundance for eligibility/ranking."""
    da=(background.control_mean_logCPM-target.control_mean_logCPM).abs()
    keep=da.le(caliper+1e-12);distance=da/caliper
    if strict:
        dd=(background[reference+'_group_dominance_share']-target[reference+'_group_dominance_share']).abs()
        keep&=background[reference+'_peak_group'].eq(target[reference+'_peak_group'])&dd.le(.2+1e-12)
        distance=distance+dd/.2
    candidates=background.loc[keep,['symbol']].copy();candidates['distance']=distance[keep]
    candidates=candidates.reset_index(drop=True).sort_values(['distance','symbol'],kind='stable')
    return candidates.symbol.tolist(),candidates.distance.tolist()


def build_covariates(root,out):
    verify_manifest(root,'replicability_snapshot');verify_manifest(root,'independent_contexts_snapshot')
    prepared=out/'prepared';prepared.mkdir(exist_ok=True);prepare(root,prepared)
    stage=pd.read_csv(root/'replicability_snapshot/all_primary_gene_stage_context.csv.gz')
    columns=['gene']+[ref+suffix for ref in REFERENCES for suffix in ['_measured','_peak_group','_peak_group_mean_CPM','_group_dominance_share']]
    stage=stage[columns].rename(columns={'gene':'symbol'})
    if stage.symbol.duplicated().any() or len(stage)!=9675:raise ValueError('Frozen universe differs')
    stats=pd.read_csv(root/'independent_contexts_snapshot/full_statistics.csv.gz')
    meta=pd.read_csv(root/'independent_contexts_snapshot/samples.csv')
    normalization=pd.read_csv(root/'independent_contexts_snapshot/normalization.csv')
    frames=[];audit=[]
    for cohort in COHORTS:
        d=stats[stats.cohort==cohort].copy()
        counts=pd.read_csv(prepared/(cohort+'.tsv.gz'),sep='\t',index_col=0)
        mm=meta[meta.cohort==cohort]
        norm=normalization[(normalization.cohort==cohort)&(normalization.variant=='full')].set_index('matrix_sample').loc[mm.matrix_sample]
        x=counts.loc[d.gene,mm.matrix_sample].to_numpy(float)
        if not np.array_equal(x.sum(axis=0),norm['lib.size'].to_numpy()):raise ValueError('Frozen filtered library sizes differ')
        effective=norm['lib.size'].to_numpy()*norm['norm.factors'].to_numpy()
        logcpm=np.log2((x+.5)/(effective+1)*1e6)
        error=float(np.max(np.abs(logcpm.mean(axis=1)-d.AveExpr.to_numpy())))
        if not np.allclose(logcpm.mean(axis=1),d.AveExpr,rtol=1e-9,atol=1e-10):raise ValueError('Voom expression reconstruction differs')
        d['control_mean_logCPM']=logcpm[:,mm.group.eq('Control').to_numpy()].mean(axis=1)
        raw_control_mean=np.log2(x[:,mm.group.eq('Control').to_numpy()]+.5).mean(axis=1)
        shift=d['control_mean_logCPM'].to_numpy()-raw_control_mean
        spread=float(np.ptp(shift))
        if spread>1e-10:raise ValueError('Unexpected gene-dependent control normalization shift')
        d=d[d.globally_unique_symbol.eq(True)&d.symbol.notna()]
        if d.symbol.duplicated().any():raise ValueError('Source symbol is not unique')
        merged=d[['cohort','gene','symbol','control_mean_logCPM','logFC','BY','BH']].merge(stage,on='symbol',how='inner',validate='one_to_one')
        frames.append(merged)
        audit.append(dict(cohort=cohort,tested_source_genes=len(stats[stats.cohort==cohort]),eligible_primary_symbols=len(merged),voom_AveExpr_max_absolute_error=error,control_libraries=int(mm.group.eq('Control').sum()),control_normalization_shift_spread=spread,normalization='Frozen full fit; source counts reconstructed from original cached assets'))
    frame=pd.concat(frames,ignore_index=True);frame.to_csv(out/'gene_covariates.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    pd.DataFrame(audit).to_csv(out/'covariate_reconstruction.csv',index=False)
    return frame


def match_contexts(covariates,priorities):
    target_rows=[];neighbor_rows=[];summary=[]
    for cohort in COHORTS:
        source=covariates[covariates.cohort==cohort].copy().set_index('symbol',drop=False)
        for reference in REFERENCES:
            valid=source.loc[stage_available(source,reference)]
            background=valid[~valid.symbol.isin(priorities)]
            for caliper in [1.0,.5]:
                matched=[]
                for symbol in sorted(priorities):
                    row=dict(cohort=cohort,reference=reference,abundance_caliper=caliper,symbol=symbol,
                             status='source_not_tested_or_nonunique',abundance_pool_size=0,stage_pool_size=0)
                    if symbol in source.index:
                        target=source.loc[symbol];row.update(control_mean_logCPM=target.control_mean_logCPM,logFC=target.logFC,
                            peak_group=target[reference+'_peak_group'],dominance_share=target[reference+'_group_dominance_share'])
                        if symbol not in valid.index:row['status']='reference_unavailable'
                        else:
                            pools={name:neighbors(target,background,reference,caliper,name!='abundance_only') for name in COMPARISONS}
                            row.update(abundance_pool_size=len(pools['abundance_only'][0]),stage_pool_size=len(pools['stage_abundance_dominance'][0]))
                            row['status']='matched' if min(row['abundance_pool_size'],row['stage_pool_size'])>=5 else 'insufficient_neighbors'
                            if row['status']=='matched':
                                matched.append(symbol)
                                for name,(ids,distances) in pools.items():
                                    for control,distance in zip(ids[:20],distances[:20]):
                                        b=background.loc[control]
                                        neighbor_rows.append(dict(cohort=cohort,reference=reference,abundance_caliper=caliper,comparison=name,
                                            target=symbol,background=control,distance=distance,
                                            abundance_gap=abs(b.control_mean_logCPM-target.control_mean_logCPM),
                                            dominance_gap=abs(b[reference+'_group_dominance_share']-target[reference+'_group_dominance_share']),
                                            target_peak_group=target[reference+'_peak_group'],background_peak_group=b[reference+'_peak_group'],
                                            target_logFC=target.logFC,background_logFC=b.logFC))
                    target_rows.append(row)
                for name in COMPARISONS:
                    summary.append(dict(cohort=cohort,reference=reference,abundance_caliper=caliper,comparison=name,
                        fixed_priorities=len(priorities),source_measured=int(source.symbol.isin(priorities).sum()),
                        reference_available=int(valid.symbol.isin(priorities).sum()),matched_priorities=len(matched),
                        annotated_background_genes=len(background)))
    return pd.DataFrame(target_rows),pd.DataFrame(neighbor_rows),pd.DataFrame(summary)


def summarize(targets,neighborhoods,summary):
    keys=['cohort','reference','abundance_caliper','comparison'];pairs=[]
    for key,d in neighborhoods.groupby(keys+['target'],sort=True):
        pairs.append(dict(zip(keys+['target'],key),target_logFC=float(d.target_logFC.iloc[0]),
            background_mean_logFC=float(d.background_logFC.mean()),background_down_1p5fold_fraction=float(d.background_logFC.le(-np.log2(1.5)).mean()),
            background_genes=len(d),maximum_abundance_gap=float(d.abundance_gap.max()),maximum_dominance_gap=float(d.dominance_gap.max())))
    per_target=pd.DataFrame(pairs,columns=keys+['target','target_logFC','background_mean_logFC','background_down_1p5fold_fraction','background_genes','maximum_abundance_gap','maximum_dominance_gap']);rows=[]
    for entry in summary.to_dict('records'):
        d=per_target
        for k in keys:d=d[d[k]==entry[k]]
        if len(d)!=entry['matched_priorities']:raise ValueError('Pair count differs')
        entry.update(priority_mean_logFC=float(d.target_logFC.mean()),background_mean_logFC=float(d.background_mean_logFC.mean()),
                     priority_down_1p5fold_fraction=float(d.target_logFC.le(-np.log2(1.5)).mean()) if len(d) else np.nan,
                     background_down_1p5fold_fraction=float(d.background_down_1p5fold_fraction.mean()))
        entry['mean_logFC_difference']=entry['priority_mean_logFC']-entry['background_mean_logFC']
        entry['down_fraction_difference']=entry['priority_down_1p5fold_fraction']-entry['background_down_1p5fold_fraction']
        excess=d.target_logFC-d.background_mean_logFC
        entry['median_pairwise_logFC_difference']=float(excess.median())
        entry['negative_pairwise_difference_fraction']=float(excess.lt(0).mean()) if len(d) else np.nan
        influence=(excess.sum()-excess)/(len(excess)-1) if len(excess)>1 else pd.Series(dtype=float)
        entry['delete_one_priority_mean_difference_min']=float(influence.min())
        entry['delete_one_priority_mean_difference_max']=float(influence.max())
        rows.append(entry)
    return per_target,pd.DataFrame(rows)


def run(root,out):
    out.mkdir(parents=True,exist_ok=True)
    for a in json.loads((root/'data/stage_specificity_inputs.json').read_text())['assets']:
        p=root/a['path']
        if p.stat().st_size!=a['bytes'] or digest(p)!=a['sha256']:raise ValueError('Input identity differs: '+a['path'])
    priorities=set(pd.read_csv(root/'replicability_snapshot/stable_descriptive_priorities.csv').gene)
    if len(priorities)!=36:raise ValueError('Priority denominator changed')
    covariates=build_covariates(root,out);targets,edges,base=match_contexts(covariates,priorities)
    pairs,summary=summarize(targets,edges,base)
    for frame,name in [(targets,'target_eligibility.csv'),(edges,'neighborhoods.csv.gz'),(pairs,'paired_comparisons.csv'),(summary,'specificity_summary.csv')]:
        frame.to_csv(out/name,index=False,compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)
    result=dict(scope='Descriptive reference-informed gene matching, not a cell-composition adjustment or inferential test',cohorts=4,reference_views=2,abundance_calipers=[1.0,.5],matching_schemes=2,summary_rows=len(summary),fixed_priorities=36,priority_eligibility_rows=len(targets),new_primary_votes=0,minimum_matchable_priorities=int(summary.matched_priorities.min()),maximum_matchable_priorities=int(summary.matched_priorities.max()),no_gene_based_Pvalues=True)
    (out/'design_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    contrasted=summary.pivot(index=['cohort','reference','abundance_caliper'],columns='comparison',values='mean_logFC_difference')
    interpretation=dict(comparison_pairs=len(contrasted),stage_mean_difference_negative=int(contrasted.stage_abundance_dominance.lt(0).sum()),
        stage_absolute_difference_smaller=int(contrasted.stage_abundance_dominance.abs().lt(contrasted.abundance_only.abs()).sum()),
        median_difference_negative=int(summary[summary.comparison=='stage_abundance_dominance'].median_pairwise_logFC_difference.lt(0).sum()),
        primary_caliper_contexts=8,no_causal_explained_fraction=True)
    (out/'interpretation_summary.json').write_text(json.dumps(interpretation,indent=2)+'\n')
    versions=dict(Python=platform.python_version(),**{p:importlib.metadata.version(p) for p in ['numpy','pandas','matplotlib']})
    (out/'runtime_versions.json').write_text(json.dumps(versions,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',default=str(ROOT/'results/stage_specificity'));args=parser.parse_args()
    run(ROOT,Path(args.output).resolve())
