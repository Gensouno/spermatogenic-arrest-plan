#!/usr/bin/env python3
"""Reconstruct original RSEM profiles; run omissions and explicitly diagnostic matched fits."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import numpy as np
import pandas as pd
from check_project import safe_path,verify_manifest

ROOT=Path(__file__).resolve().parents[1]
REQUIRED=['full','exclude_pair']+[f'omit_GSM{i}' for i in range(7034750,7034790)]


def check_inputs(root):
    for a in json.loads((root/'data/setx_stability_inputs.json').read_text())['assets']:
        p=safe_path(root,a['path'])
        if p.stat().st_size!=a['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=a['sha256']:
            raise ValueError('Stability input identity differs: '+a['path'])


def prepare(root,out):
    check_inputs(root)
    for name in ['setx_interaction_snapshot','setx_length_snapshot']:verify_manifest(root,name)
    audit=pd.read_csv(root/'new_models_snapshot/Setx_archive_audit.csv').sort_values('gsm')
    columns={k:{} for k in ['expected_count','effective_length','TPM']};genes=None
    with tarfile.open(root/'data/raw/replicability_discovery/GSE224875_RAW.tar') as tar:
        for row in audit.itertuples():
            payload=tar.extractfile(row.member).read()
            if hashlib.sha256(payload).hexdigest()!=row.member_sha256:raise ValueError('Genome member hash differs')
            frame=pd.read_csv(io.BytesIO(gzip.decompress(payload)),sep='\t')
            if genes is None:genes=frame.gene_id.copy()
            if genes.isna().any() or genes.duplicated().any() or not genes.equals(frame.gene_id):raise ValueError('Gene IDs differ')
            for key in columns:
                values=frame[key].to_numpy(float)
                if not np.isfinite(values).all() or (values<0).any():raise ValueError('Invalid RSEM field')
                columns[key][row.gsm]=values
    for key,values in columns.items():
        pd.DataFrame(values,index=pd.Index(genes,name='gene')).to_csv(out/(key+'.tsv.gz'),sep='\t',compression={'method':'gzip','mtime':0})


def candidate_stability(effects,candidates,priorities,required=REQUIRED):
    if effects.duplicated(['variant','symbol']).any():raise ValueError('Duplicate candidate variant/symbol')
    if set(effects.variant.unique())!=set(required):raise ValueError('Incorrect or incomplete variant suite')
    rows=[]
    for symbol in candidates:
        sub=effects[effects.symbol==symbol].set_index('variant')
        complete=set(required)<=set(sub.index)
        rows.append(dict(symbol=symbol,fixed_stable_priority=symbol in priorities,variants_measured=len(sub),
            variants_required=len(required),all_variants_measured=complete,
            minimum_delta_P18=float(sub.delta_P18.min()) if len(sub) else np.nan,
            maximum_delta_P18=float(sub.delta_P18.max()) if len(sub) else np.nan,
            maximum_delta_variant=str(sub.delta_P18.idxmax()) if len(sub) else '',
            all_negative=bool(complete and (sub.delta_P18<0).all()),
            all_negative_1p5fold=bool(complete and (sub.delta_P18<=-np.log2(1.5)).all()),
            conditional_BY05_all_variants_descriptive=bool(complete and (sub.BY<.05).all()),
            main_BY=float(sub.loc['full','BY']) if 'full' in sub.index else np.nan))
    return pd.DataFrame(rows)


def summarize(root,out):
    effects=pd.read_csv(out/'candidate_omissions.csv.gz')
    candidates=pd.read_csv(root/'replicability_snapshot/candidate_stage_context.csv').gene.tolist()
    priorities=set(pd.read_csv(root/'replicability_snapshot/stable_descriptive_priorities.csv').gene)
    stable=candidate_stability(effects,candidates,priorities)
    prior=pd.read_csv(root/'setx_interaction_snapshot/candidate_stability.csv')
    stable=stable.merge(prior[['symbol','all_variants_down_1p5fold']].rename(columns={'all_variants_down_1p5fold':'prior_TPM_all_variants_1p5fold'}),on='symbol',how='left',validate='one_to_one')
    stable['both_count_and_TPM_suites_1p5fold']=stable.all_negative_1p5fold&stable.prior_TPM_all_variants_1p5fold.fillna(False)
    stable.to_csv(out/'candidate_stability.csv',index=False)
    diagnostics=pd.read_csv(out/'matched_diagnostics.csv.gz');summaries=[]
    for name,d in diagnostics.groupby('variant',sort=True):
        priority=d[d.symbol.isin(priorities)]
        summaries.append(dict(variant=name,tested_genes=len(d),conditional_BY05=int((d.BY<.05).sum()),conditional_BH05=int((d.BH<.05).sum()),
            fixed_priorities_measured=len(priority),fixed_priorities_negative=int((priority.delta_P18<0).sum()),
            fixed_priorities_negative_1p5fold=int((priority.delta_P18<=-np.log2(1.5)).sum()),
            fixed_priorities_conditional_BY05=int((priority.BY<.05).sum())))
    pd.DataFrame(summaries).to_csv(out/'matched_summary.csv',index=False)
    sets={name:set(d.loc[d.symbol.isin(priorities)&(d.BY<.05),'symbol']) for name,d in diagnostics.groupby('variant')}
    length_set=sets['matched_count_length_offset'];zero_set=sets['matched_count_zero_offset_counterfactual']
    transitions=dict(fixed_priorities=len(priorities),length_offset_only=sorted(length_set-zero_set),
        zero_offset_counterfactual_only=sorted(zero_set-length_set),both=sorted(length_set&zero_set),
        scope='Threshold membership on the fixed common universe; no additional subgroup FDR guarantee or preferred zero-offset workflow.')
    (out/'matched_threshold_transitions.json').write_text(json.dumps(transitions,indent=2)+'\n')
    p=stable[stable.fixed_stable_priority]
    report=dict(scope='Post-hoc same-library stability and diagnostic comparison; no biological confirmation.',
        fixed_candidates=len(candidates),fixed_priorities=len(priorities),count_variants=len(REQUIRED),single_library_omissions=40,
        priorities_complete=int(p.all_variants_measured.sum()),priorities_always_negative=int(p.all_negative.sum()),
        priorities_always_negative_1p5fold=int(p.all_negative_1p5fold.sum()),
        priorities_conditional_BY05_all_variants_descriptive=int(p.conditional_BY05_all_variants_descriptive.sum()),
        priorities_1p5fold_both_count_and_TPM_suites=int(p.both_count_and_TPM_suites_1p5fold.sum()),
        matched_models=len(summaries),matched_tested_genes=int(summaries[0]['tested_genes']),
        caveats=['A missing/re-filtered/non-unique candidate fails retention.','Omission ranges are not confidence intervals or independent replication.',
            'Per-variant adjusted-threshold intersections and candidate subsets have no additional FDR guarantee.',
            'Zero-length-offset counterfactual is a diagnostic, not a recommended estimate-count workflow.',
            'Matched-count fits share normalization but reestimate dispersions; this does not fully attribute the original TPM/count discrepancy.',
            'Animal independence, pooling, source code 7844 and barcode nesting remain unresolved; no primary genetic-unit vote added.'])
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    plot_stability(root,out,stable,diagnostics,priorities)
    print(json.dumps(report,indent=2),flush=True)
    print(pd.read_csv(out/'matched_summary.csv').to_string(index=False),flush=True)


def plot_stability(root,out,stable,diagnostics,priorities):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    order=sorted(priorities);rows=stable.set_index('symbol').reindex(order)
    effects=pd.read_csv(out/'candidate_omissions.csv.gz')
    baseline=effects[effects.variant=='full'].set_index('symbol').reindex(order)
    tpm=pd.read_csv(root/'setx_interaction_snapshot/candidate_stability.csv').set_index('symbol').reindex(order)
    fig,axes=plt.subplots(1,2,figsize=(11,11),gridspec_kw={'width_ratios':[2.5,1]})
    for i,symbol in enumerate(order):
        row=rows.loc[symbol]
        if np.isfinite(row.maximum_delta_P18):
            axes[0].plot([row.minimum_delta_P18,row.maximum_delta_P18],[i,i],color='.55',lw=1)
            axes[0].plot(baseline.loc[symbol,'delta_P18'],i,'o',color='tab:blue',ms=3,
                label='Count full fit' if i==0 else None)
            axes[0].plot(tpm.loc[symbol,'main_delta_P18'],i,'o',mfc='none',mec='tab:red',ms=4,
                label='Prior TPM full fit' if i==0 else None)
    axes[0].set_yticks(range(len(order)),order,fontsize=7)
    axes[0].set_ylim(len(order)-.5,-.5);axes[0].axvline(0,color='.3',lw=.7)
    axes[0].axvline(-np.log2(1.5),color='tab:red',ls=':',lw=.8)
    axes[0].set_title('P18-minus-P10 change in genotype contrast',fontsize=10)
    axes[0].set_xlabel('Count range over 42 fits; range is not a confidence interval',fontsize=8)
    axes[0].legend(loc='upper left',fontsize=7)
    names=['matched_TPM_limma','matched_count_length_offset','matched_count_zero_offset_counterfactual']
    columns=[]
    for name in names:
        q=diagnostics[(diagnostics.variant==name)&diagnostics.symbol.isin(order)].set_index('symbol').reindex(order).BY
        columns.append(np.where(q.notna(),(q<.05).astype(float),np.nan))
    values=np.column_stack(columns);cmap=ListedColormap(['#eeeeee','#2166ac']);cmap.set_bad('#999999')
    axes[1].imshow(values,aspect='auto',cmap=cmap,vmin=0,vmax=1)
    axes[1].set_yticks([]);axes[1].set_xticks(range(3),['TPM\nlimma','Count\nlength offset','Count\nzero offset'],fontsize=8,rotation=45,ha='right')
    axes[1].set_title('Conditional BY < 0.05\nSame 16,853-gene universe',fontsize=10)
    axes[1].set_xlabel('Blue: met; light gray: not met; dark gray: missing\nZero offset is a diagnostic counterfactual',fontsize=8)
    fig.suptitle('Setx: fixed 36 descriptive priorities\nSame-library exploratory sensitivities; biological independence unresolved',fontsize=11)
    fig.subplots_adjust(left=.15,right=.97,top=.91,bottom=.14,wspace=.3)
    fig.savefig(out/'stability_and_matched_models.png',dpi=160);plt.close(fig)


def main():
    out=ROOT/'results/setx_stability';out.mkdir(exist_ok=True,parents=True)
    prepare(ROOT,out)
    subprocess.run(['bash','envs/rscript.sh','code/analyze_setx_stability.R'],cwd=ROOT,check=True)
    summarize(ROOT,out)


if __name__=='__main__':main()
