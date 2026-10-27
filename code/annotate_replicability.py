#!/usr/bin/env python3
"""Independent expression context for every primary gene, never a discovery filter."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from check_stage_reference import read_centroids,FILENAME
from analyze_replicability import PRIMARY,UNITS

ROOT=Path(__file__).resolve().parents[1]
CHEN={'Spermatogonia':['A1','In','BS','BG2'],
      'Preleptotene':['G1','ePL','mPL','lPL'],
      'Meiotic':['L','Z','eP','mP','lP','D','MI','MII'],
      'Postmeiotic':['RS2','RS4','RS6','RS8']}
GREEN={'Spermatogonia':['GC1'],'Meiotic_or_transition':[f'GC{i}' for i in range(2,9)],
       'Postmeiotic':[f'GC{i}' for i in range(9,13)]}


def summarize_context(profile,genes,groups,prefix):
    """Mean stages equally; retain missing, zeros and tied maxima explicitly."""
    if profile.index.has_duplicates or not np.isfinite(profile.to_numpy()).all() or (profile<0).any().any():
        raise ValueError('Invalid source profile')
    columns=[c for stages in groups.values() for c in stages]
    if len(columns)!=len(set(columns)) or set(columns)!=set(profile.columns):raise ValueError('Incomplete or overlapping stage groups')
    source=profile.reindex(genes)
    broad=pd.DataFrame({g:source[stages].mean(axis=1) for g,stages in groups.items()})
    measured=source.notna().all(axis=1);expressed=measured & (source.sum(axis=1)>0)
    peak=broad.max(axis=1);share=peak/broad.sum(axis=1).replace(0,np.nan)
    tied=np.isclose(broad.to_numpy(),peak.to_numpy()[:,None],rtol=1e-8,atol=1e-10).sum(axis=1)>1
    finepeak=source.max(axis=1)
    finetied=np.isclose(source.to_numpy(),finepeak.to_numpy()[:,None],rtol=1e-8,atol=1e-10).sum(axis=1)>1
    label=pd.Series('not_measured',index=source.index)
    fine_label=label.copy()
    # Only assign an argmax to measured, expressed, non-tied rows.
    good=expressed & ~tied
    label.loc[good]=broad.loc[good].idxmax(axis=1)
    label.loc[expressed & tied]='tied_groups'
    label.loc[measured & ~expressed]='unexpressed_in_source'
    goodfine=expressed & ~finetied
    fine_label.loc[goodfine]=source.loc[goodfine].idxmax(axis=1)
    fine_label.loc[expressed & finetied]='tied_stages'
    fine_label.loc[measured & ~expressed]='unexpressed_in_source'
    result=pd.DataFrame({prefix+'_measured':measured,prefix+'_peak_group':label,
        prefix+'_peak_stage':fine_label,prefix+'_peak_group_mean_CPM':peak,
        prefix+'_group_dominance_share':share,
        prefix+'_high_dominance_descriptive':expressed & ~tied & (peak>=10) & (share>=.5)})
    for g in groups:result[prefix+'_'+g+'_mean_CPM']=broad[g]
    return result


def main():
    out=ROOT/'results/replicability'
    gene=pd.read_csv(out/'primary_gene_summary.csv.gz',index_col='gene')
    transfer=json.loads((ROOT/'reference_limits_snapshot/transfer_design.json').read_text())
    fixed=transfer['fixed_prepared_inputs_sha256']
    paths=['results/reference_limits/chen_stage_profiles.csv.gz','results/reference_limits/chen_profile_evidence.csv',
           'results/reference_limits/chen_preparation_audit.json']
    fingerprints={}
    for relative in paths:
        h=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        if h!=fixed[relative]:raise ValueError('Frozen Chen source differs: '+relative)
        fingerprints[relative]=h
    relative='data/raw/ref/'+FILENAME
    expected=json.loads((ROOT/'results_snapshot/stage_reference/provenance.json').read_text())['sha256']
    h=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
    if h!=expected:raise ValueError('Frozen Green workbook differs')
    fingerprints[relative]=h
    chen=pd.read_csv(ROOT/paths[0],index_col='gene')
    chen=chen[[f'Chen:{s}:equal_cell_mean' for stages in CHEN.values() for s in stages]]
    chen.columns=[c.split(':')[1] for c in chen]
    green,cells,identity_error,excluded=read_centroids(ROOT/relative)
    context=summarize_context(chen,gene.index,CHEN,'Chen').join(summarize_context(green,gene.index,GREEN,'Green'))
    comparable={'Spermatogonia':'Spermatogonia','Meiotic':'Meiotic_or_transition','Postmeiotic':'Postmeiotic'}
    context['coarse_source_agreement']='not_comparable'
    eligible=context.Chen_peak_group.isin(comparable) & context.Green_peak_group.isin(GREEN)
    context.loc[eligible,'coarse_source_agreement']=np.where(
        context.loc[eligible,'Chen_peak_group'].map(comparable)==context.loc[eligible,'Green_peak_group'],'agree','disagree')
    # Add all five cohort effects, retaining both Rfx2 ages as separate columns.
    study=pd.read_csv(out/'recomputed_study_statistics.csv.gz')
    effects=study[study.gse.isin(PRIMARY)].pivot(index='gene',columns='gse',values='logFC').reindex(gene.index)
    context=context.join(effects.add_suffix('_logFC'))
    context.to_csv(out/'all_primary_gene_stage_context.csv.gz',compression={'method':'gzip','mtime':0})
    candidate=gene[(gene.unit_up_1p5fold>=3)|(gene.unit_down_1p5fold>=3)].join(context)
    candidate['descriptive_direction']=np.where(candidate.unit_up_1p5fold>=3,'up','down')
    candidate['LOO_stable_at_least_3_units']=np.where(candidate.descriptive_direction=='up',
        candidate.LOO_stable_up_1p5fold_units>=3,candidate.LOO_stable_down_1p5fold_units>=3)
    unit_loo=pd.read_csv(out/'LOO_gene_unit_summary.csv.gz').set_index(['gene','unit'])
    support=[];baseline_names=[];stable_names=[]
    for symbol,row in candidate.iterrows():
        direction=row.descriptive_direction
        base=[];retained=[]
        for unit in UNITS:
            effect=row[unit+('_logFC_min' if direction=='up' else '_logFC_max')]*(1 if direction=='up' else -1)
            fits=effect>=np.log2(1.5)
            omission=unit_loo.loc[(symbol,unit)]
            stable=bool(omission['stable_'+direction+'_1p5fold'])
            if fits:base.append(unit)
            if fits and stable:retained.append(unit)
            support.append({'gene':symbol,'direction':direction,'unit':unit,'directional_minimum_log2FC':effect,
                'baseline_1p5fold_support':fits,'all_LOO_variants_available':bool(omission.all_variants_available),
                'LOO_1p5fold_support':stable,'baseline_and_LOO_support':fits and stable,
                'directional_P':row[unit+'_P_'+direction]})
        baseline_names.append('/'.join(base));stable_names.append('/'.join(retained))
    candidate['baseline_supporting_units'] = baseline_names
    candidate['LOO_supporting_units'] = stable_names
    candidate['baseline_and_LOO_stable_units']=[len(n.split('/')) if n else 0 for n in stable_names]
    candidate['LOO_stable_at_least_3_units']=candidate.baseline_and_LOO_stable_units>=3
    pd.DataFrame(support).to_csv(out/'candidate_unit_support.csv',index=False)
    magnitudes=[]
    for _,row in candidate.iterrows():
        v=[row[u+('_logFC_min' if row.descriptive_direction=='up' else '_logFC_max')]
           *(1 if row.descriptive_direction=='up' else -1) for u in UNITS]
        magnitudes.append(sorted(v,reverse=True)[2])
    candidate['third_largest_directional_unit_log2FC']=magnitudes
    candidate['status']='Descriptive follow-up priority only; no replicability discovery or validated mechanism'
    candidate.to_csv(out/'candidate_stage_context.csv')
    stable=candidate[candidate.LOO_stable_at_least_3_units].copy()
    stable=stable.sort_values(['Chen_peak_group','third_largest_directional_unit_log2FC'],ascending=[True,False])
    stable.to_csv(out/'stable_descriptive_priorities.csv')
    rows=[];cross=[]
    for subset,frame in [('all_primary',context),('descriptive_candidates',candidate),('LOO_stable_candidates',stable)]:
        for (a,b),n in frame.groupby(['Chen_peak_group','Green_peak_group']).size().items():
            cross.append({'subset':subset,'Chen_peak_group':a,'Green_peak_group':b,'genes':int(n),'denominator':len(frame)})
        for source in ['Chen','Green']:
            for label,n in frame[source+'_peak_group'].value_counts().items():
                rows.append({'subset':subset,'source':source,'peak_group':label,'genes':int(n),'denominator':len(frame)})
        for label,n in frame.coarse_source_agreement.value_counts().items():
            rows.append({'subset':subset,'source':'coarse_source_agreement','peak_group':label,'genes':int(n),'denominator':len(frame)})
    pd.DataFrame(rows).to_csv(out/'stage_context_summary.csv',index=False)
    pd.DataFrame(cross).to_csv(out/'stage_source_cross_counts.csv',index=False)
    design={'source_sha256':fingerprints,'Chen_groups':CHEN,'Green_groups':GREEN,
        'source_genes':{'Chen':len(chen),'Green_valid_symbols':len(green),'Green_excluded_nonstring_rows':len(excluded)},
        'Green_natural_log_identity_error':identity_error,
        'primary_genes_annotated':len(context),'descriptive_candidates':len(candidate),'LOO_stable_candidates':len(stable),
        'aggregation':'Chen equal-cell-normalized CPM; Green linearized mean normalized CPM. Equal-stage/centroid arithmetic means within each predefined group, no renormalization over selected genes.',
        'matching':'Exact symbols only; missing, measured all-zero, tied groups/stages and incomparable preleptotene explicit.',
        'high_dominance_flag':'Unique peak, group mean >=10 CPM and peak/sum(group means)>=0.5. Descriptive only; no validated marker or inferential filtering.',
        'limitations':['Stage samples/centroids are expression context, not independent animals for a gene-enrichment test.',
            'Both sources are germ-cell views; missing somatic types, age/protocol and symbol differences limit interpretation.',
            'Green GC2-GC8 is meiotic-or-transition, not exact leptotene/zygotene.',
            'No cell-number proportions, causal composition mediation or cell-intrinsic shared mechanism inferred.']}
    (out/'stage_annotation_design.json').write_text(json.dumps(design,indent=2)+'\n')
    print(json.dumps(design,indent=2),flush=True)


if __name__=='__main__':main()
