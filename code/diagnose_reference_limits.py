#!/usr/bin/env python3
"""Prespecified gene-universe, weighting and generic bulk-length sensitivities."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

from check_stage_reference import read_centroids,collapse_centroids
from diagnose_reference_transfer import MODELS,broad_weights,gene_folds,train_fit
from diagnose_leptotene import identifier_map,read_htseq,split_htseq_rows,COMPARATORS

ROOT=Path(__file__).resolve().parents[1]


def merged_length(intervals):
    if not intervals:return 0
    pairs=sorted(intervals)
    if any(a<1 or b<a for a,b in pairs):raise ValueError('Invalid one-based inclusive interval')
    total=0;start,stop=pairs[0]
    for a,b in pairs[1:]:
        if a<=stop+1:stop=max(stop,b)
        else:total+=stop-start+1;start,stop=a,b
    return total+stop-start+1


def exon_lengths(path,mapping):
    intervals={};locations={};rows=0
    with gzip.open(path,'rt') as stream:
        for line in stream:
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)!=9:raise ValueError('Invalid GTF field count')
            if fields[2]!='exon':continue
            match=re.search(r'(?:^|;\s*)gene_id "(ENSMUSG[0-9]+)"',fields[8])
            if not match or match[1] not in mapping.index:continue
            gene=match[1];where=(fields[0],fields[6]);rows+=1
            if gene in locations and locations[gene]!=where:raise ValueError('Mapped gene has multiple chromosome/strand locations')
            locations[gene]=where;intervals.setdefault(gene,[]).append((int(fields[3]),int(fields[4])))
    lengths=pd.Series({mapping.loc[g]:merged_length(v) for g,v in intervals.items()},name='union_exon_bp',dtype=float)
    if lengths.index.duplicated().any() or (lengths<=0).any():raise ValueError('Ambiguous/nonpositive gene length')
    return lengths,{'source_exon_rows':rows,'unique_gene_lengths':len(lengths),'definition':'Union of all annotated exons per gene, one-based inclusive, no double-counting transcript overlap.'}


def extended_folds(legacy,expanded):
    old=gene_folds(sorted(legacy));new=sorted(set(expanded)-set(old))
    old.update(dict(zip(new,np.random.default_rng(20261025).permutation(len(new))%5)))
    return old


def normalize(x):
    if x.index.has_duplicates or not np.isfinite(x.to_numpy()).all() or (x<0).any().any() or (x.sum()<=0).any():raise ValueError('Invalid abundance matrix')
    return x.div(x.sum(),axis=1)*1e6


def full_references():
    path=ROOT/'data/raw/ref/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx'
    audit=json.loads((ROOT/'results/stage_reference/provenance.json').read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest()!=audit['sha256']:raise ValueError('Original Green centroid workbook changed')
    fine,cells,error,excluded=read_centroids(path)
    full_valid_genes=len(fine)
    pooled=pd.read_csv(ROOT/'results/reference_calibration/pooled_cpm.csv.gz',index_col=0).T
    # Independently reconstruct the pooled profile from the old cached class sums.
    from calibrate_reference_batches import signature,CLASSES
    data=np.load(ROOT/'data/ref/batches/counts.npz');meta=pd.read_csv(ROOT/'data/ref/batches/batch_classes.csv')
    rebuilt=pd.DataFrame(signature(data['sums'],meta,np.ones(len(meta),dtype=bool)),index=data['genes'],columns=CLASSES)
    pd.testing.assert_frame_equal(pooled,rebuilt,check_names=False,check_index_type=False,check_column_type=False,rtol=1e-10,atol=1e-9)
    genes=fine.index.intersection(pooled.index).sort_values();fine=fine.loc[genes];broad=collapse_centroids(fine,cells)
    detailed=fine.copy();detailed.insert(0,'Somatic',pooled.loc[genes,'Somatic'])
    broad.insert(0,'Somatic',pooled.loc[genes,'Somatic'])
    refs={'original_4':pooled.loc[genes],'centroid_4':broad,'centroid_13':detailed}
    old={m:pd.read_csv(ROOT/f'results/stage_reference/pooled_{m}_reference.csv.gz',index_col=0) for m in MODELS}
    for m in MODELS:
        test=normalize(refs[m].loc[old[m].index]).reindex(columns=old[m].columns)
        pd.testing.assert_frame_equal(test,old[m],check_names=False,check_index_type=False,check_column_type=False,rtol=1e-10,atol=1e-8)
    return refs,old,{'full_green_valid_genes':full_valid_genes,'full_Green_pooled_common_genes':len(genes),'legacy_reference_genes':len(old['centroid_4']),
        'excel_nonstring_gene_rows_excluded':len(excluded),'original_Average_Total_identity_max_error':error,
        'legacy_reference_reconstruction_matches':True,'pooled_reference_reconstructed_from_cached_class_sums':True}


def load_targets():
    vectors=[];rows=[];mapping=identifier_map(pd.read_csv(ROOT/'data/raw/labeled_reference/genes.tsv',sep='\t'))
    spec=json.loads((ROOT/'data/leptotene_assets.json').read_text())
    evidence=pd.read_csv(ROOT/'leptotene_snapshot/sample_evidence.csv').set_index('profile')
    for r in spec['selections']:
        vector,_=read_htseq(ROOT/r['path'],mapping);name=r['gsm'];vectors.append(vector.rename(name))
        rows.append({'profile':name,'target_set':'GSE155967','author_type':r['stage'],'target_weighting':'RNA_counts',
            'expected_broad_class':evidence.loc[name,'expected_broad_class']})
    mixed,_=split_htseq_rows(pd.read_csv(ROOT/'data/raw/leptotene/expression/GSE83264_mRNA_genes_merged_counts.txt.gz',sep='\t',index_col=0))
    for name,(stage,col) in COMPARATORS.items():
        vectors.append(mixed[col].rename(name));rows.append({'profile':name,'target_set':'GSE83264','author_type':stage,
            'target_weighting':'RNA_counts','expected_broad_class':'Meiotic'})
    x=pd.read_csv(ROOT/'results/emptydrops_p15/type_profiles.csv.gz',index_col=0)
    e=pd.read_csv(ROOT/'emptydrops_p15_snapshot/profile_evidence.csv')
    e=e[(e.branch=='EmptyDrops')&(e.stratum=='full')]
    if len(e)!=30:raise ValueError('Expected all 15 full-P15 groups and both weighting views')
    for r in e.itertuples():
        vectors.append(x[r.profile]);rows.append({'profile':r.profile,'target_set':'E-MTAB-6946_P15','author_type':r.author_type,
            'target_weighting':r.weighting,'expected_broad_class':r.expected_broad_class})
    chen=ROOT/'results/reference_limits/chen_stage_profiles.csv.gz'
    if chen.exists():
        x=pd.read_csv(chen,index_col=0);e=pd.read_csv(chen.with_name('chen_profile_evidence.csv'))
        if x.columns.tolist()!=e.profile.tolist():raise ValueError('Chen profile identities differ')
        for r in e.itertuples():
            vectors.append(x[r.profile]);rows.append({'profile':r.profile,'target_set':'GSE107644','author_type':r.author_type,
                'target_weighting':r.weighting,'expected_broad_class':r.expected_broad_class})
    return vectors,pd.DataFrame(rows),mapping


def main():
    out=ROOT/'results/reference_limits';out.mkdir(parents=True,exist_ok=True)
    refs,old,source_audit=full_references();vectors,evidence,mapping=load_targets()
    # No outcome-based feature selection; all admitted targets share these sets.
    common=refs['centroid_4'].index
    for vector in vectors:common=common.intersection(vector.index)
    common=common.sort_values();legacy=common.intersection(old['centroid_4'].index).sort_values()
    if len(common)<=len(legacy) or len(legacy)<10000:raise ValueError('Unexpected gene-universe coverage')
    gtf=ROOT/'data/ref/annotation/Mus_musculus.GRCm38.79.gtf.gz'
    expected=json.loads((ROOT/'data/ref/annotation/provenance.json').read_text())['sha256']
    if hashlib.sha256(gtf.read_bytes()).hexdigest()!=expected:raise ValueError('Frozen legacy GTF differs')
    lengths,length_audit=exon_lengths(gtf,mapping);lengths.to_csv(out/'gene_union_exon_lengths.csv.gz',index_label='gene',compression={'method':'gzip','mtime':0})
    folds=extended_folds(old['centroid_4'].index,refs['centroid_4'].index)
    metrics=[];weights=[];conditions=[];gene_rows=[];prediction_rows=[]
    for universe,genes in [('legacy',legacy),('expanded',common)]:
        for transform in ['raw','bulk_length_covered_raw','bulk_union_exon_sensitivity']:
            use=genes if transform=='raw' else genes.intersection(lengths.index).sort_values()
            selected=evidence if transform=='raw' else evidence[evidence.target_set=='GSE155967']
            source=[v for v in vectors if v.name in set(selected.profile)]
            target=pd.concat([v.loc[use] for v in source],axis=1)
            if transform=='bulk_union_exon_sensitivity':target=target.div(lengths.loc[use],axis=0)
            y=normalize(target);references={m:normalize(r.loc[use]) for m,r in refs.items()}
            fold=np.array([folds[g] for g in use])
            for weighting in ['inverse_sqrt','uniform']:
                scale=1/np.sqrt(references['centroid_4'].mean(axis=1).to_numpy()+1) if weighting=='inverse_sqrt' else np.ones(len(use))
                design={'gene_universe':universe,'target_transform':transform,'gene_weighting':weighting,'genes':len(use)}
                gene_rows.extend({'gene':g,'fold':int(f),'weight':float(w),**design} for g,f,w in zip(use,fold,scale))
                for model,reference in references.items():
                    s=reference.to_numpy();conditions.append({**design,'model':model,'weighted_condition_number':float(np.linalg.cond(s*scale[:,None]))})
                    for row in selected.to_dict('records'):
                        t=y[row['profile']].to_numpy();prediction=np.empty_like(t);scored=pd.notna(row['expected_broad_class'])
                        for f in [-1,0,1,2,3,4]:
                            coef=train_fit(s,t,scale,np.ones(len(use),dtype=bool) if f==-1 else fold!=f)
                            if f!=-1:prediction[fold==f]=s[fold==f]@coef
                            p=broad_weights(coef,reference.columns,model)
                            weights.append({**design,**row,'model':model,'fold':f,'largest_class':p.idxmax(),
                                'matches_author_broad':p.idxmax()==row['expected_broad_class'] if scored else None,
                                'expected_class_weight':p[row['expected_broad_class']] if scored else None,**p.to_dict()})
                        original=next(v for v in source if v.name==row['profile'])
                        metrics.append({**design,**row,'model':model,'source_expression_mass_on_genes':float(original.loc[use].sum()/original.sum()),
                            'holdout_log2_RMSE':float(np.sqrt(np.mean((np.log2(t+.5)-np.log2(prediction+.5))**2))),
                            'holdout_weighted_RMSE':float(np.sqrt(np.mean(((t-prediction)*scale)**2)))})
                print('Completed',universe,transform,weighting,len(use),'genes;',len(selected),'targets',flush=True)
    w=pd.DataFrame(weights);m=pd.DataFrame(metrics);full=w[w.fold==-1]
    w.to_csv(out/'technical_fold_weights.csv.gz',index=False,compression={'method':'gzip','mtime':0});m.to_csv(out/'technical_holdout_metrics.csv',index=False)
    full.to_csv(out/'technical_full_fit_weights.csv',index=False)
    keys=['profile','model','gene_universe','target_transform','gene_weighting']
    full.merge(m[keys+['holdout_log2_RMSE','holdout_weighted_RMSE','source_expression_mass_on_genes']],on=keys,validate='one_to_one').to_csv(out/'technical_comparison.csv',index=False)
    pd.DataFrame(conditions).to_csv(out/'technical_conditions.csv',index=False);pd.DataFrame(gene_rows).to_csv(out/'technical_gene_design.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    evidence.to_csv(out/'target_evidence.csv',index=False)
    source_paths=['data/raw/ref/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx','results/reference_calibration/pooled_cpm.csv.gz',
        'data/ref/batches/counts.npz','data/ref/batches/batch_classes.csv','results/stage_reference/provenance.json',
        'data/ref/annotation/Mus_musculus.GRCm38.79.gtf.gz','data/ref/annotation/provenance.json','data/raw/labeled_reference/genes.tsv',
        'results/emptydrops_p15/type_profiles.csv.gz','emptydrops_p15_snapshot/profile_evidence.csv','leptotene_snapshot/sample_evidence.csv']
    source_paths+=[r['path'] for r in json.loads((ROOT/'data/leptotene_assets.json').read_text())['assets']]
    source_paths += [f'results/stage_reference/pooled_{m}_reference.csv.gz' for m in MODELS]
    if (out/'chen_stage_profiles.csv.gz').exists():source_paths+=['results/reference_limits/chen_stage_profiles.csv.gz','results/reference_limits/chen_profile_evidence.csv']
    design={'scope':'Gene-universe, gene-weighting and generic union-exon bulk sensitivity; no target-driven marker selection or calibrated cell fractions.',
        'source_reference_audit':source_audit,'target_profiles':len(evidence),'expanded_common_genes':len(common),'legacy_common_genes':len(legacy),
        'length_annotation':{'release':'Ensembl 79 GRCm38','source_url':'https://ftp.ensembl.org/pub/release-79/gtf/mus_musculus/Mus_musculus.GRCm38.79.gtf.gz',**length_audit,
            'limitation':'Generic same-build historical annotation, not asserted the original GTF or validated effective transcript length. No gene-length division of UMI targets or UMI references; no mm10 lengths applied to mm9 comparator.'},
        'folds':'Legacy seed 20261024 assignments retained; new sorted reference genes appended with seed 20261025, before target intersection.',
        'normalization':'Each common-gene reference/target column normalized before gene folds. Different universes and transformations change target values and weights; cross-universe errors are not on an identical response scale. Within-universe model comparisons share the same response.',
        'length_comparison':'bulk_length_covered_raw and bulk_union_exon_sensitivity use identical genes, reference normalization, weights and folds; only the target count/length transformation changes. These are expression diagnostics, not original TPM reconstruction.',
        'selection':'Exploratory diagnosis motivated by previous L/Z findings; not an untouched independent validation.',
        'fixed_inputs_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(set(source_paths))}}
    (out/'technical_design.json').write_text(json.dumps(design,indent=2)+'\n')


if __name__=='__main__':main()
