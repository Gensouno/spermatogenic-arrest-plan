#!/usr/bin/env python3
"""All fixed candidates in new expression contexts, with explicit missingness."""
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from check_replicability_context import unique_effects,transfer_summary
ROOT=Path(__file__).resolve().parents[1]


def mapped_effects(frame,mapping):
    if frame.gene.duplicated().any() or mapping.gene.duplicated().any():raise ValueError('Ambiguous source gene identifiers')
    d=frame.merge(mapping,on='gene',how='left',validate='one_to_one')
    d.loc[d.symbol.isin(['-','--','']), 'symbol']=np.nan
    return unique_effects(d,'symbol')


def main():
    out=ROOT/'results/new_models';base=ROOT/'results/replicability'
    candidates=pd.read_csv(base/'candidate_stage_context.csv')
    primary=pd.read_csv(base/'primary_gene_summary.csv.gz').gene
    fitted=pd.read_csv(out/'conditional_model_statistics.csv.gz')
    loo=pd.read_csv(out/'leave_one_out_effects.csv.gz')
    pseudo=pd.read_csv(out/'TPM_pseudocount_effects.csv.gz')
    annotation=pd.read_csv(out/'gene_annotation.csv')
    rows=[];summaries=[];coverage=[]
    cohorts=['Eef1g_LZ_P12']+[f'Setx_P{i}' for i in [10,12,14,16,18]]
    if set(fitted.cohort)!=set(cohorts):raise ValueError('New contrast inventory differs')
    for cohort in cohorts:
        source='Eef1g' if cohort.startswith('Eef1g') else 'Setx'
        mapping=annotation[annotation.cohort_source==source][['gene','symbol']]
        de=fitted[fitted.cohort==cohort]
        effect=mapped_effects(de,mapping)
        omitted=loo[loo.cohort==cohort]
        omissions={g:mapped_effects(p,mapping) for g,p in omitted.groupby('deleted')}
        if len(omissions)!=(6 if source=='Eef1g' else 8):raise ValueError('Omission count differs')
        part=transfer_summary(effect,omissions,candidates.gene,np.where(candidates.descriptive_direction=='up',1,-1))
        part['cohort']=cohort;part['unit']='counts' if source=='Eef1g' else 'TPM'
        part['primary_candidate_direction']=candidates.descriptive_direction.to_numpy()
        part['primary_LOO_priority']=candidates.LOO_stable_at_least_3_units.to_numpy()
        part['pseudocount_sensitivity_applicable']=source=='Setx'
        if source=='Setx':
            alternatives={v:mapped_effects(p,mapping).reindex(candidates.gene).to_numpy()
                for v,p in pseudo[pseudo.cohort==cohort].groupby('pseudocount')}
            if set(alternatives)!={.1,1}:raise ValueError('Pseudocount inventory differs')
            signed=np.stack([part.candidate_directional_external_log2FC,
                *[alternatives[v]*np.where(candidates.descriptive_direction=='up',1,-1) for v in [.1,1]]],axis=1)
            part['all_pseudocount_1p5fold']=np.isfinite(signed).all(axis=1)&(signed>=np.log2(1.5)).all(axis=1)
            part['all_LOO_and_pseudocount_1p5fold']=part.baseline_and_all_LOO_1p5fold & part.all_pseudocount_1p5fold
        else:
            part['all_pseudocount_1p5fold']=pd.NA;part['all_LOO_and_pseudocount_1p5fold']=pd.NA
        rows.append(part)
        for subset,frame in [('all_105_candidates',part),('primary_36_priorities',part[part.primary_LOO_priority])]:
            for direction in ['all','up','down']:
                d=frame if direction=='all' else frame[frame.primary_candidate_direction==direction]
                summaries.append({'cohort':cohort,'unit':part.unit.iloc[0],'subset':subset,'direction':direction,'fixed_genes':len(d),
                    **{k:int(d[k].sum()) for k in ['tested_exact_unique_symbol','sign_matches','directional_1p5fold',
                        'all_LOO_available','baseline_and_all_LOO_sign_match','baseline_and_all_LOO_1p5fold']},
                    'all_LOO_and_pseudocount_1p5fold':int(d.all_LOO_and_pseudocount_1p5fold.sum()) if source=='Setx' else None})
        coverage.append({'cohort':cohort,'tested_gene_features':len(de),'tested_unique_symbols':len(effect),
            'primary_genes':len(primary),'primary_covered':int(effect.reindex(primary).notna().sum()),
            'candidate_genes':len(candidates),'candidate_covered':int(part.tested_exact_unique_symbol.sum()),
            'not_unique_or_unmapped_tested_features':len(de)-len(effect)})
    combined=pd.concat(rows,ignore_index=True)
    combined.to_csv(out/'candidate_context.csv',index=False)
    summary=pd.DataFrame(summaries);summary.to_csv(out/'context_summary.csv',index=False)
    pd.DataFrame(coverage).to_csv(out/'gene_coverage.csv',index=False)
    fixed=['results/replicability/candidate_stage_context.csv','results/replicability/primary_gene_summary.csv.gz']
    design={'source_contract':'data/new_model_audit_contract.json','fixed_prior_inputs_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in fixed},
        'contrasts':cohorts,'primary_genetic_votes_added':0,'fixed_candidates':len(candidates),'fixed_priorities':int(candidates.LOO_stable_at_least_3_units.sum()),
        'filter':'Eef1g filterByExpr before TMM/voom/limma; Setx TPM>=1 in at least the smaller group size, log2(TPM+0.5), limma trend and robust. Re-filter in every omission.',
        'annotation':'Eef1g original provided gene_name; Setx current installed org.Mm.eg.db unique ENSEMBL->SYMBOL. Fit every retained original gene ID; exclude all duplicated symbols separately within each full or omitted fit for cross-study mapping.',
        'multiplicity':'Conditional model adjusted P-values apply to source gene features, not confirmatory candidate-list replication. No biological FDR or cell-intrinsic mechanism claim.',
        'limits':['Public processed expression, not raw-read remapping or source DESeq2 reproduction.',
            'Replicate titles do not establish independent animals, pooling or hidden batches; model P-values are conditional on those assumptions.',
            'Setx five ages are one gene model; do not infer age-specific interactions from separate fits.',
            'Eef1g FACS L/Z mixture is not whole testis or pure fixed stage; sorting and state makeup remain potential confounders.',
            'No causal composition explanation or validation of prior deconvolution.']}
    (out/'design.json').write_text(json.dumps(design,indent=2)+'\n')
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    stable=candidates[candidates.LOO_stable_at_least_3_units].sort_values('gene')
    matrix=combined.pivot(index='gene',columns='cohort',values='external_logFC').reindex(index=stable.gene,columns=cohorts)
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd')
    fig,ax=plt.subplots(figsize=(7,10),layout='constrained');im=ax.imshow(matrix,aspect='auto',cmap=cmap,vmin=-2,vmax=2)
    ax.set_yticks(range(len(matrix)),matrix.index,fontsize=8)
    ax.set_xticks(range(6),['Eef1g L/Z\nP12 counts',*[f'Setx P{i}\nTPM' for i in [10,12,14,16,18]]],rotation=50,ha='right',fontsize=8)
    ax.set_title('Fixed 36 downward priorities in new contexts\nKO - source control; grey = not tested / unique mapped',fontsize=10)
    fig.colorbar(im,ax=ax,shrink=.5,label='log2FC; colour clipped at +/-2')
    fig.savefig(out/'fixed_priorities_context.png',dpi=180);plt.close(fig)
    pca=pd.read_csv(out/'PCA.csv');fig,axes=plt.subplots(2,3,figsize=(11,6),layout='constrained')
    for ax,cohort in zip(axes.flat,cohorts):
        d=pca[pca.cohort==cohort]
        for group,part in d.groupby('group'):
            ax.scatter(part.PC1,part.PC2,label=group,s=30)
            for i,(_,r) in enumerate(part.iterrows()):
                ax.annotate(r.gsm[-2:],(r.PC1,r.PC2),fontsize=6,xytext=(4,5 if i%2 else -8),textcoords='offset points')
        ax.set(title=cohort,xlabel=f'PC1 {d.PC1_variance.iloc[0]:.1%}',ylabel=f'PC2 {d.PC2_variance.iloc[0]:.1%}')
        ax.legend(frameon=False,fontsize=7)
        ax.margins(x=.18,y=.16)
    fig.suptitle('Exploratory sample QC; source libraries, no automatic exclusions')
    fig.savefig(out/'sample_QC.png',dpi=180);plt.close(fig)
    print(summary[(summary.subset=='primary_36_priorities')&(summary.direction=='all')].to_string(index=False),flush=True)


if __name__=='__main__':main()
