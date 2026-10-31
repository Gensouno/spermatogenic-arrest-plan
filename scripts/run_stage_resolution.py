#!/usr/bin/env python3
"""Post-hoc joint matching of native group-mean and individual-stage peaks."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from run_stage_specificity import ROOT,COHORTS,REFERENCES,stage_available,neighbors,summarize

SCHEMES=['abundance_only','stage_abundance_dominance','stage_peak_abundance_dominance']


def descriptor_audit(stage,priorities,design):
    rows=[]
    for ref in REFERENCES:
        mapping={label:group for group,labels in design[ref+'_groups'].items() for label in labels}
        for subset,d in [('all_primary',stage),('fixed_priorities',stage[stage.symbol.isin(priorities)])]:
            valid=stage_available(d,ref)&d[ref+'_peak_stage'].isin(mapping)
            v=d[valid];implied=v[ref+'_peak_stage'].map(mapping)
            rows.append(dict(reference=ref,subset=subset,total=len(d),comparable=len(v),
                peak_stage_in_other_mean_peak_group=int(implied.ne(v[ref+'_peak_group']).sum()),
                meaning='Different aggregation descriptors can legitimately disagree; no relabeling or data correction.'))
    return pd.DataFrame(rows)


def run(root,out):
    cov=pd.read_csv(out/'gene_covariates.csv.gz')
    stage=pd.read_csv(root/'replicability_snapshot/all_primary_gene_stage_context.csv.gz').rename(columns={'gene':'symbol'})
    priorities=set(pd.read_csv(root/'replicability_snapshot/stable_descriptive_priorities.csv').gene)
    stagecols=['symbol','Chen_peak_stage','Green_peak_stage']
    cov=cov.merge(stage[stagecols],on='symbol',how='left',validate='many_to_one')
    descriptor_audit(stage,priorities,json.loads((root/'replicability_snapshot/stage_annotation_design.json').read_text())).to_csv(out/'stage_descriptor_audit.csv',index=False)
    targets=[];edges=[];base=[]
    for cohort in COHORTS:
        source=cov[cov.cohort==cohort].set_index('symbol',drop=False)
        for ref in REFERENCES:
            valid=source.loc[stage_available(source,ref)&source[ref+'_peak_stage'].notna()]
            bg=valid[~valid.symbol.isin(priorities)]
            for caliper in [1.,.5]:
                matched=[]
                for symbol in sorted(priorities):
                    row=dict(cohort=cohort,reference=ref,abundance_caliper=caliper,symbol=symbol,status='source_or_reference_unavailable',fine_pool_size=0)
                    if symbol in valid.index:
                        t=valid.loc[symbol]
                        pools={s:neighbors(t,bg if s!='stage_peak_abundance_dominance' else bg[bg[ref+'_peak_stage'].eq(t[ref+'_peak_stage'])],ref,caliper,s!='abundance_only') for s in SCHEMES}
                        row.update(fine_pool_size=len(pools['stage_peak_abundance_dominance'][0]),status='matched' if min(len(p[0]) for p in pools.values())>=5 else 'insufficient_neighbors')
                        if row['status']=='matched':
                            matched.append(symbol)
                            for scheme,(ids,distances) in pools.items():
                                for b,distance in zip(ids[:20],distances[:20]):
                                    q=bg.loc[b]
                                    edges.append(dict(cohort=cohort,reference=ref,abundance_caliper=caliper,comparison=scheme,target=symbol,background=b,distance=distance,
                                        abundance_gap=abs(q.control_mean_logCPM-t.control_mean_logCPM),dominance_gap=abs(q[ref+'_group_dominance_share']-t[ref+'_group_dominance_share']),
                                        target_peak_group=t[ref+'_peak_group'],background_peak_group=q[ref+'_peak_group'],
                                        target_peak_stage=t[ref+'_peak_stage'],background_peak_stage=q[ref+'_peak_stage'],target_logFC=t.logFC,background_logFC=q.logFC))
                    targets.append(row)
                for scheme in SCHEMES:base.append(dict(cohort=cohort,reference=ref,abundance_caliper=caliper,comparison=scheme,fixed_priorities=36,
                    source_measured=int(source.symbol.isin(priorities).sum()),reference_available=int(valid.symbol.isin(priorities).sum()),matched_priorities=len(matched),annotated_background_genes=len(bg)))
    eligibility=pd.DataFrame(targets);edge=pd.DataFrame(edges);pairs,summary=summarize(eligibility,edge,pd.DataFrame(base))
    eligibility.to_csv(out/'resolution_eligibility.csv',index=False);edge.to_csv(out/'resolution_neighborhoods.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    pairs.to_csv(out/'resolution_pairs.csv',index=False);summary.to_csv(out/'resolution_summary.csv',index=False)
    wide=summary.pivot(index=['cohort','reference','abundance_caliper'],columns='comparison',values='mean_logFC_difference')
    report=dict(scope='Post-hoc native peak-stage constraint on the same targets within each scenario, not independent validation.',summary_rows=len(summary),comparison_pairs=len(wide),
        minimum_matchable_priorities=int(summary.matched_priorities.min()),maximum_matchable_priorities=int(summary.matched_priorities.max()),
        fine_negative_mean=int(wide.stage_peak_abundance_dominance.lt(0).sum()),
        fine_absolute_difference_smaller_than_coarse=int(wide.stage_peak_abundance_dominance.abs().lt(wide.stage_abundance_dominance.abs()).sum()),
        fine_negative_median=int(summary[summary.comparison=='stage_peak_abundance_dominance'].median_pairwise_logFC_difference.lt(0).sum()),
        no_new_primary_votes=True)
    (out/'resolution_design.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default=str(ROOT/'results/stage_specificity'));a=p.parse_args();run(ROOT,Path(a.output))
