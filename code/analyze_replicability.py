#!/usr/bin/env python3
"""Four genetic-unit directional replication, with explicit descriptive stability."""
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from replicability_statistics import checked_pvalues,intersection_union,partial_conjunction,adjust_joint_family,directional_tail_from_two_sided

ROOT=Path(__file__).resolve().parents[1]
PRIMARY=['GSE236388','GSE283904','GSE68283','GSE68284','GSE101511']
SENSITIVITY=['GSE294209','GSE38096','GSE42601']
UNITS={'XAP5':['GSE236388'],'N6amt1':['GSE283904'],'Rfx2':['GSE68283','GSE68284'],'Hnrnpgt':['GSE101511']}


def tested_common(table,cohorts):
    return sorted(set.intersection(*(set(table.loc[table.gse==g,'gene']) for g in cohorts)))


def make_unit_evidence(table,genes):
    pieces={g:table[table.gse==g].set_index('gene').loc[genes] for g in PRIMARY}
    frame=pd.DataFrame(index=pd.Index(genes,name='gene'));up=[];down=[];pos=[];neg=[]
    for unit,cohorts in UNITS.items():
        p_up=intersection_union(np.stack([pieces[g].P_up for g in cohorts],axis=1))
        p_down=intersection_union(np.stack([pieces[g].P_down for g in cohorts],axis=1))
        effects=np.stack([pieces[g].logFC for g in cohorts],axis=1)
        low=effects.min(axis=1);high=effects.max(axis=1)
        frame[unit+'_P_up']=p_up;frame[unit+'_P_down']=p_down
        frame[unit+'_logFC_min']=low;frame[unit+'_logFC_max']=high
        frame[unit+'_direction']=np.where(low>0,'up',np.where(high<0,'down',np.where((low==0)&(high==0),'zero','mixed_or_zero')))
        up.append(p_up);down.append(p_down);pos.append(low);neg.append(high)
    return frame,{'up':np.stack(up,axis=1),'down':np.stack(down,axis=1)},np.stack(pos,axis=1),np.stack(neg,axis=1)


def run_universe(table,genes,label):
    model,pvalues,pos,neg=make_unit_evidence(table,genes)
    records=[]
    for direction in ['up','down']:
        for r in [3,4]:
            p=partial_conjunction(pvalues[direction],r)
            records.extend({'gene':g,'universe':label,'direction':direction,'r':r,'genetic_units':4,
                'partial_conjunction_P':float(value)} for g,value in zip(genes,p))
    tests=pd.DataFrame(records)
    # A single joint family, never separately adjust favourable directions/levels.
    tests['joint_BY_q']=adjust_joint_family(tests.partial_conjunction_P.to_numpy(),'by')
    tests['joint_BH_q_sensitivity']=adjust_joint_family(tests.partial_conjunction_P.to_numpy(),'bh')
    tests['BY_rejected']=tests.joint_BY_q<=.05;tests['BH_rejected_sensitivity']=tests.joint_BH_q_sensitivity<=.05
    model['unit_up_count']=(pos>0).sum(axis=1);model['unit_down_count']=(neg<0).sum(axis=1)
    for threshold,name in [(np.log2(1.5),'1p5fold'),(1.0,'2fold')]:
        model['unit_up_'+name]=(pos>=threshold).sum(axis=1)
        model['unit_down_'+name]=(neg<=-threshold).sum(axis=1)
    return tests,model


def loo_summaries(root,genes):
    table=pd.read_csv(root/'results/replicability/recomputed_count_LOO_effects.csv.gz')
    table=table[table.gse.isin(PRIMARY)].copy()
    if table.duplicated(['gse','variant','gene']).any() or not np.isfinite(table.logFC).all():raise ValueError('Invalid cached leave-one-out table')
    expected={g:len(pd.read_csv(root/f'data/bulk/{g}_samples.tsv',sep='\t')) for g in PRIMARY}
    observed=table.groupby('gse').variant.nunique().to_dict()
    if observed!=expected or sum(expected.values())!=40:raise ValueError('Leave-one-out variant/sample reconciliation differs')
    table=table[table.gene.isin(genes)]
    table['positive']=table.logFC>0;table['negative']=table.logFC<0
    summary=table.groupby(['gene','gse']).agg(available_variants=('variant','size'),
        min_logFC=('logFC','min'),max_logFC=('logFC','max'),positive_variants=('positive','sum'),negative_variants=('negative','sum'))
    summary=summary.reindex(pd.MultiIndex.from_product([genes,PRIMARY],names=['gene','gse'])).reset_index()
    summary['expected_variants']=summary.gse.map(expected)
    summary['available_variants']=summary.available_variants.fillna(0).astype(int)
    summary['all_variants_available']=summary.available_variants==summary.expected_variants
    unit_rows=[]
    for unit,cohorts in UNITS.items():
        for gene,part in summary[summary.gse.isin(cohorts)].groupby('gene',sort=False):
            complete=bool(part.all_variants_available.all())
            low=float(part.min_logFC.min()) if part.available_variants.sum()>0 else np.nan
            high=float(part.max_logFC.max()) if part.available_variants.sum()>0 else np.nan
            unit_rows.append({'gene':gene,'unit':unit,'all_variants_available':complete,
                'available_cohort_variants':int(part.available_variants.sum()),'expected_cohort_variants':int(part.expected_variants.sum()),
                'min_logFC':low,'max_logFC':high,'stable_up':complete and low>0,'stable_down':complete and high<0,
                'stable_up_1p5fold':complete and low>=np.log2(1.5),'stable_down_1p5fold':complete and high<=-np.log2(1.5)})
    units=pd.DataFrame(unit_rows)
    counts=units.groupby('gene')[['stable_up','stable_down','stable_up_1p5fold','stable_down_1p5fold']].sum()
    counts.columns=['LOO_'+c+'_units' for c in counts]
    return summary,units,counts


def main():
    out=ROOT/'results/replicability';out.mkdir(parents=True,exist_ok=True)
    (out/'runtime_Python_versions.json').write_text(json.dumps({'Python':platform.python_version(),
        'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},indent=2)+'\n')
    table=pd.read_csv(out/'recomputed_study_statistics.csv.gz')
    if table.duplicated(['gse','gene']).any() or table.gene.isna().any() or set(table.gse)!=set(PRIMARY+SENSITIVITY):raise ValueError('Source statistic inventory differs')
    for c in ['P_up','P_down','P_two_sided','original_adj_P']:checked_pvalues(table[c])
    if not np.isfinite(table[['logFC','t','df_total']]).all().all() or (table.df_total<=0).any():raise ValueError('Invalid fitted statistics')
    for direction,column in [('up','P_up'),('down','P_down')]:
        np.testing.assert_allclose(directional_tail_from_two_sided(table.P_two_sided,table.t,direction),table[column],atol=1e-12,rtol=1e-10)
    primary=tested_common(table,PRIMARY);legacy=tested_common(table,PRIMARY+SENSITIVITY)
    previous=pd.read_csv(ROOT/'results/all_cohort_comparison/common_logFC.csv.gz',index_col=0)
    if len(legacy)!=8099 or set(legacy)!=set(previous.index):raise ValueError('Legacy gene-universe reconciliation differs')
    for gse in PRIMARY+SENSITIVITY:
        np.testing.assert_allclose(previous[gse].reindex(legacy),table[table.gse==gse].set_index('gene').logFC.reindex(legacy),atol=1e-8,rtol=1e-8)
    alltests=[];models={}
    for name,genes in [('counts_primary',primary),('legacy_eight_cohort_sensitivity',legacy)]:
        tests,model=run_universe(table,genes,name);alltests.append(tests);models[name]=model
        model.to_csv(out/(name+'_unit_evidence.csv.gz'),compression={'method':'gzip','mtime':0})
    tests=pd.concat(alltests,ignore_index=True);tests.to_csv(out/'directional_replicability_tests.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    primary_tests=tests[tests.universe=='counts_primary'].copy();models['counts_primary'].to_csv(out/'primary_gene_summary.csv.gz',compression={'method':'gzip','mtime':0})
    findings=primary_tests[primary_tests.BY_rejected].copy()
    findings=findings.merge(models['counts_primary'].reset_index(),on='gene',validate='many_to_one')
    findings.to_csv(out/'primary_BY_findings.csv',index=False)
    primary_tests.nsmallest(100,'partial_conjunction_P').assign(status='Ranked for inspection; use BY_rejected field, not rank as discovery').to_csv(out/'top_ranked_hypotheses.csv',index=False)
    loo,unit_loo,stable=loo_summaries(ROOT,primary)
    loo.to_csv(out/'LOO_gene_cohort_summary.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    unit_loo.to_csv(out/'LOO_gene_unit_summary.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    gene=models['counts_primary'].join(stable)
    # Cross-platform effects are annotations only; never primary votes.
    for gse in SENSITIVITY:
        other=table[table.gse==gse].set_index('gene')
        gene[gse+'_logFC_sensitivity']=other.logFC.reindex(gene.index)
    qwide=primary_tests.pivot(index='gene',columns=['direction','r'],values='joint_BY_q')
    for direction in ['up','down']:
        for r in [3,4]:gene[f'BY_q_{direction}_{r}of4_joint_family']=qwide[(direction,r)]
    gene.to_csv(out/'primary_gene_summary.csv.gz',compression={'method':'gzip','mtime':0})
    descriptive=gene[(gene.unit_up_1p5fold>=3)|(gene.unit_down_1p5fold>=3)].copy()
    descriptive['interpretation']='Effect-size/direction screen only; not a new FDR-controlled list or validated mechanism'
    descriptive.to_csv(out/'descriptive_1p5fold_candidates.csv.gz',compression={'method':'gzip','mtime':0})
    pattern=gene[[u+'_direction' for u in UNITS]].agg('/'.join,axis=1)
    pattern.value_counts().rename_axis('XAP5_N6amt1_Rfx2_Hnrnpgt').rename('genes').to_csv(out/'direction_patterns.csv')
    summary=tests.groupby(['universe','direction','r']).agg(tested_hypotheses=('gene','size'),
        BY_rejections=('BY_rejected','sum'),BH_rejections_sensitivity=('BH_rejected_sensitivity','sum'),
        minimum_raw_P=('partial_conjunction_P','min'),minimum_BY_q=('joint_BY_q','min')).reset_index()
    summary.to_csv(out/'replicability_summary.csv',index=False)
    desc=[]
    for universe,model in models.items():
        for name,threshold in [('sign',0),('1p5fold',np.log2(1.5)),('2fold',1)]:
            for direction in ['up','down']:
                key='unit_'+direction+'_count' if threshold==0 else 'unit_'+direction+'_'+name
                for r in [3,4]:desc.append({'universe':universe,'direction':direction,'r':r,'effect_threshold':name,
                    'genes_meeting_direction_effect_screen':int((model[key]>=r).sum()),'tested_genes':len(model),
                    'inferential_status':'Descriptive support only, no FDR guarantee for this selected subset'})
    pd.DataFrame(desc).to_csv(out/'descriptive_support_counts.csv',index=False)
    coverage=[]
    for gse,part in table.groupby('gse'):
        coverage.append({'gse':gse,'unit':part.unit.iloc[0],'tested_genes':len(part),'in_primary_common':int(part.gene.isin(primary).sum()),
            'in_legacy_common':int(part.gene.isin(legacy).sum()),'role':'Primary count cohort' if gse in PRIMARY else 'Sensitivity annotation only'})
    pd.DataFrame(coverage).to_csv(out/'gene_coverage.csv',index=False)
    design={'scope':'Four genetic-unit directional expression replication; whole-testis effects only.',
        'primary_count_cohorts':PRIMARY,'genetic_units':UNITS,'count_samples':40,'primary_common_tested_genes':len(primary),'legacy_common_tested_genes':len(legacy),
        'primary_hypothesis_family_size':len(primary_tests),'directions':['up','down'],'replicability_levels':[3,4],
        'Rfx2':'Maximum of age-specific directional P-values; both ages required for that unit. Sign/effect screens use age-wise min for positive and max for negative, never an average hiding disagreement.',
        'test':'Bonferroni partial conjunction min(1,(n-r+1)*p_(r)), n=4; valid under arbitrary component dependence if component P-values are valid.',
        'multiplicity':'BY across all genes x both directions x both r levels jointly within each universe. BH is a sensitivity requiring additional dependence assumptions. Primary claims refer to hypothesis-level family, not a separately guaranteed FDR for a filtered/deduplicated gene subset.',
        'Pvalue_source':'Moderated-t one-sided pt() tails from recomputed source fits; original logFC/t/two-sided/adjusted P-values reconciled first.',
        'effect_screen':'Thresholds 0, log2(1.5), 1 are separate descriptive summaries. Do not use an effect or LOO filter to imply a new FDR guarantee.',
        'LOO':'All 40 count one-sample-deletion variants recomputed from prepared matrices and reconciled with the prior cache; filtering changes can remove genes. Completeness required before marking a unit stable. Not independent replication or uncertainty intervals.',
        'limitations':['Source limma/filtering assumptions and P-value validity are required; no empirical biological FDR calibration.',
            'Limited Hnrnpgt power and unconfirmed N6amt1 age remain.',
            'Bulk effects can reflect cellular makeup; no cell-intrinsic mechanism, composition-adjusted effect or severity ranking inferred.',
            'No new raw preprocessing, transcript annotation reconstruction or array normalization in this branch.'],
        'primary_BY_discovery_hypotheses':len(findings),'primary_BY_unique_genes':int(findings.gene.nunique()),
        'descriptive_1p5fold_candidate_genes':len(descriptive)}
    (out/'design.json').write_text(json.dumps(design,indent=2)+'\n')
    print(summary.to_string(index=False),flush=True);print(json.dumps(design,indent=2),flush=True)


if __name__=='__main__':main()
