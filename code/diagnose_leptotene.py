#!/usr/bin/env python3
"""Independent stage-separated and mixed-LZ profiles against unchanged references."""
import gzip
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from diagnose_reference_transfer import MODELS,broad_weights,gene_folds,train_fit
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from discover_expansion import parse_soft
from fetch_leptotene import STAGES

ROOT=Path(__file__).resolve().parents[1]
COMPARATORS={'GSM2197870':('LZ','LZn1WT'),'GSM2197871':('LZ','LZn2WT'),
    'GSM2197872':('Pachytene','Pachn1WT'),'GSM2197873':('Pachytene','Pachn2WT')}


def checked_counts(table):
    values=table.to_numpy(dtype=float)
    if table.index.has_duplicates or table.index.isna().any() or not np.isfinite(values).all() or (values<0).any() or not (values==np.floor(values)).all():
        raise ValueError('Expected unique gene identifiers and nonnegative integer counts')
    if (values.sum(axis=0)<=0).any():raise ValueError('Empty count profile')
    return table


def identifier_map(genes):
    if genes.ID.isna().any() or genes.ID.duplicated().any():raise ValueError('Ambiguous ID map')
    mouse=genes.ID.str.fullmatch(r'ENSMUSG[0-9]+')
    keep=mouse & genes.Symbol.notna() & ~genes.loc[mouse,'Symbol'].duplicated(False).reindex(genes.index,fill_value=False)
    return genes.loc[keep].set_index('ID').Symbol


def split_htseq_rows(table):
    checked_counts(table)
    summary=table.index.str.startswith('__')
    known={'__no_feature','__ambiguous','__too_low_aQual','__not_aligned','__alignment_not_unique'}
    if not set(table.index[summary]).issubset(known):raise ValueError('Unexpected HTSeq summary key')
    genes=table.loc[~summary]
    checked_counts(genes)
    return genes,table.loc[summary]


def read_htseq(path,mapping):
    table=pd.read_csv(path,sep='\t',header=None,index_col=0)
    if table.shape[1]!=1:raise ValueError('Expected two-column HTSeq count file')
    genes,summary=split_htseq_rows(table)
    if not genes.index.str.fullmatch(r'ENSMUSG[0-9]+').all():raise ValueError('Expected unversioned mouse Ensembl gene IDs')
    keep=genes.index.isin(mapping.index);mapped=genes.loc[keep].iloc[:,0].copy();mapped.index=mapping.loc[mapped.index].to_numpy()
    if mapped.index.duplicated().any():raise ValueError('Mapped symbol ambiguity')
    audit={'gene_rows':len(genes),'mapped_gene_rows':int(keep.sum()),'unmapped_or_ambiguous_gene_rows':int((~keep).sum()),
        'gene_count_total':int(genes.to_numpy().sum()),'mapped_count_total':int(mapped.sum()),
        'mapping_expression_mass':float(mapped.sum()/genes.to_numpy().sum()),
        'excluded_HTSeq_summary_counts':{str(k):int(v) for k,v in summary.iloc[:,0].items()}}
    return mapped,audit


def metadata(gse):
    path=ROOT/f'data/raw/leptotene/metadata/{gse}_family.soft.gz'
    series,samples=parse_soft(gzip.decompress(path.read_bytes()).decode())
    if len(series.get('Series_sample_id',[]))!=len(samples):raise ValueError('Metadata sample count differs')
    return series,samples


def paper_evidence(path):
    root=ET.parse(path).getroot()
    dois=[e.text for e in root.iter('article-id') if e.get('pub-id-type')=='doi']
    if '10.1038/s41467-021-26033-0' not in dois:raise ValueError('Wrong primary paper XML')
    sections={s.findtext('title',''):' '.join(s.itertext()) for s in root.iter('sec')}
    text=sections.get('RNA-seq library construction','')
    if not all(word in text for word in ['10-day-old','16-day-old','2-week-old','two independent biological samples']):
        raise ValueError('Stage age/replicate evidence not found in primary source')
    return {'doi':'10.1038/s41467-021-26033-0','source_url':'https://www.nature.com/articles/s41467-021-26033-0',
        'Preleptotene':'10-day-old unsynchronized mice','Leptotene':'16-day-old synchronized mice, WIN 18,446 then RA',
        'Zygotene':'2-week-old unsynchronized mice','replication':'Paper describes two independent biological samples per RNA stage; animal pooling and exact per-sample IDs not resolved.',
        'genotype':'Primary methods describe C57BL/6J and a separate SUN1 W151R experiment; RNA sample metadata do not contain an explicit genotype field. Do not promote to genotype-verified WT validation.',
        'purity':'Paper reports enriched sorted populations, not 100% pure cell-fraction truth.'}


def main():
    out=ROOT/'results/leptotene_followup';out.mkdir(parents=True,exist_ok=True)
    spec=json.loads((ROOT/'data/leptotene_assets.json').read_text())
    for row in spec['assets']+[json.loads((ROOT/'data/leptotene_paper_asset.json').read_text())]:
        path=ROOT/row['path']
        if path.is_symlink() or path.stat().st_size!=row['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
            raise ValueError('Frozen raw input differs')
    mapping_path=ROOT/'data/raw/labeled_reference/genes.tsv'
    original=next(r for r in json.loads((ROOT/'data/labeled_reference_assets.json').read_text())['assets'] if Path(r['path']).name=='genes.tsv')
    if hashlib.sha256(mapping_path.read_bytes()).hexdigest()!=original['sha256']:raise ValueError('Frozen gene mapping source differs')
    mapping=identifier_map(pd.read_csv(mapping_path,sep='\t'))
    paper=ROOT/'data/raw/leptotene/paper/PMC8501046.xml';provenance=paper_evidence(paper)
    evidence=[];vectors=[];mapping_audits=[];_,samples=metadata('GSE155967')
    if {k for k,v in samples.items() if v.get('Sample_library_strategy')==['RNA-Seq']}!=set(STAGES):raise ValueError('RNA inventory changed')
    selected={s['gsm']:s for s in spec['selections']}
    for gsm,stage in STAGES.items():
        record=samples[gsm]
        if record['Sample_organism_ch1']!=['Mus musculus'] or 'cell type: '+stage not in record['Sample_characteristics_ch1']:
            raise ValueError('Stage/species mismatch')
        filename=Path(selected[gsm]['path']).name
        declared=[u for k,v in record.items() if k.startswith('Sample_supplementary_file') for u in v]
        if len(declared)!=1 or Path(declared[0]).name!=filename:raise ValueError('RNA source identity mismatch')
        vector,audit=read_htseq(ROOT/selected[gsm]['path'],mapping);vectors.append(vector.rename(gsm))
        mapping_audits.append({'profile':gsm,'gse':'GSE155967',**audit})
        evidence.append({'profile':gsm,'gse':'GSE155967','stage':stage,'expected_broad_class':'Meiotic' if stage!='Preleptotene' else None,
            'unit':'HTSeq gene counts','genome_build':'mm10','original_title':record['Sample_title'][0],
            'preparation':provenance[stage],'genotype_status':'not explicitly specified per RNA sample',
            'original_column':filename,'mapping_expression_mass':audit['mapping_expression_mass']})
    _,samples=metadata('GSE83264');matrix,summary=split_htseq_rows(pd.read_csv(ROOT/'data/raw/leptotene/expression/GSE83264_mRNA_genes_merged_counts.txt.gz',sep='\t',index_col=0))
    for gsm,(stage,column) in COMPARATORS.items():
        r=samples[gsm];rep='n1' if column.endswith('n1WT') else 'n2'
        cell_type='purified leptotene/zygotene spermatocytes' if stage=='LZ' else 'purified pachytene spermatocytes'
        expected=['genotype/variation: wild type','molecule subtype: RNA','replicate: '+rep,'tissue/cell type: '+cell_type]
        if r['Sample_organism_ch1']!=['Mus musculus'] or r.get('Sample_library_strategy')!=['RNA-Seq'] or not all(s in r['Sample_characteristics_ch1'] for s in expected):raise ValueError('Mixed-LZ comparator identity mismatch')
        if not any('70-80' in text for text in r.get('Sample_growth_protocol_ch1',[])):raise ValueError('Adult comparator age evidence differs')
        vector=matrix[column];vectors.append(vector.rename(gsm))
        mapping_audits.append({'profile':gsm,'gse':'GSE83264','gene_rows':len(matrix),'mapped_gene_rows':len(matrix),
            'unmapped_or_ambiguous_gene_rows':0,'gene_count_total':int(vector.sum()),'mapped_count_total':int(vector.sum()),
            'mapping_expression_mass':1.0,'excluded_HTSeq_summary_counts':{str(k):int(v) for k,v in summary[column].items()}})
        evidence.append({'profile':gsm,'gse':'GSE83264','stage':stage,'expected_broad_class':'Meiotic','unit':'HTSeq gene counts',
            'genome_build':'mm9','original_title':r['Sample_title'][0],'preparation':'Purified adult 70-80 day spermatocytes; LZ is a mixed population',
            'genotype_status':'explicit wild type','original_column':column,'mapping_expression_mass':1.0})
    # Common to every admitted source and model, never selected by fit quality.
    genes=vectors[0].index
    for vector in vectors[1:]:genes=genes.intersection(vector.index)
    refpaths=[f'results/stage_reference/pooled_{m}_reference.csv.gz' for m in MODELS]
    refs={m:pd.read_csv(ROOT/p,index_col=0) for m,p in zip(MODELS,refpaths)}
    prior=json.loads((ROOT/'labeled_reference_snapshot/design.json').read_text())['fixed_inputs_sha256']
    if any(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=prior[p] for p in refpaths):raise ValueError('Reference modified')
    genes=genes.intersection(refs['centroid_4'].index).sort_values()
    if len(genes)<10000:raise ValueError('Unexpectedly low common gene coverage')
    evidence=pd.DataFrame(evidence)
    for vector in vectors:
        evidence.loc[evidence.profile.eq(vector.name),'common_gene_expression_mass_of_mapped']=vector.loc[genes].sum()/vector.sum()
    x=pd.concat([v.loc[genes] for v in vectors],axis=1);x.index.name='gene'
    x.to_csv(out/'common_gene_counts.csv.gz');evidence.to_csv(out/'sample_evidence.csv',index=False)
    pd.DataFrame(mapping_audits).to_json(out/'mapping_audit.json',orient='records',indent=2)
    y=x.div(x.sum(),axis=1)*1e6
    refs={m:r.loc[genes].div(r.loc[genes].sum(),axis=1)*1e6 for m,r in refs.items()}
    foldmap=gene_folds(sorted(pd.read_csv(ROOT/refpaths[1],index_col=0).index))
    folds=np.array([foldmap[g] for g in genes]);scale=1/np.sqrt(refs['centroid_4'].mean(axis=1).to_numpy()+1)
    weights=[];metrics=[];conditions=[]
    for model,ref in refs.items():
        s=ref.to_numpy();conditions.append({'model':model,'genes':len(genes),'weighted_condition_number':float(np.linalg.cond(s*scale[:,None]))})
        for row in evidence.to_dict('records'):
            target=y[row['profile']].to_numpy();prediction=np.empty_like(target);scored=pd.notna(row['expected_broad_class'])
            for fold in [-1,0,1,2,3,4]:
                coef=train_fit(s,target,scale,np.ones(len(genes),dtype=bool) if fold==-1 else folds!=fold)
                if fold!=-1:prediction[folds==fold]=s[folds==fold]@coef
                p=broad_weights(coef,ref.columns,model)
                weights.append({**row,'model':model,'fold':fold,'largest_class':p.idxmax(),
                    'matches_author_broad':p.idxmax()==row['expected_broad_class'] if scored else None,
                    'expected_class_weight':p[row['expected_broad_class']] if scored else None,**p.to_dict()})
            metrics.append({**row,'model':model,'genes':len(genes),
                'holdout_log2_RMSE':float(np.sqrt(np.mean((np.log2(target+.5)-np.log2(prediction+.5))**2))),
                'holdout_weighted_RMSE':float(np.sqrt(np.mean(((target-prediction)*scale)**2)))})
        print('Diagnosed',model,len(evidence),'profiles;',len(genes),'fixed common genes',flush=True)
    weights=pd.DataFrame(weights);metrics=pd.DataFrame(metrics);full=weights[weights.fold==-1]
    weights.to_csv(out/'fold_weights.csv',index=False);full.to_csv(out/'full_fit_weights.csv',index=False);metrics.to_csv(out/'holdout_metrics.csv',index=False)
    full.merge(metrics[['profile','model','holdout_log2_RMSE','holdout_weighted_RMSE']],on=['profile','model'],validate='one_to_one').to_csv(out/'stage_comparison.csv',index=False)
    pd.DataFrame(conditions).to_csv(out/'design_conditions.csv',index=False)
    pd.DataFrame({'gene':genes,'fold':folds,'weight':scale}).to_csv(out/'gene_design.csv.gz',index=False)
    replicates=[]
    for (gse,stage),rows in evidence.groupby(['gse','stage']):
        if len(rows)!=2:raise ValueError('Expected two reported samples per stage/source')
        a,b=rows.profile.tolist();replicates.append({'gse':gse,'stage':stage,'first':a,'second':b,
            'common_genes':len(genes),'Spearman_log2_CPM':float(np.log2(y[a]+.5).corr(np.log2(y[b]+.5),method='spearman'))})
    pd.DataFrame(replicates).to_csv(out/'replicate_concordance.csv',index=False)
    barcode=pd.read_csv(ROOT/'data/raw/leptotene/expression/GSE107644_barcode_information.txt.gz',sep='\t')
    review=[{'gse':'GSE155967','samples_in_GEO':22,'admitted_RNA_profiles':6,'assessable_broad_labels':4,'decision':'Stage-separated external diagnosis; preL unscored, synchronization and incomplete per-sample genotype evidence retained.'},
        {'gse':'GSE83264','samples_in_GEO':20,'admitted_RNA_profiles':4,'assessable_broad_labels':4,'decision':'Independent mixed LZ and pachytene comparators only; six whole-testis RNA and ten small-RNA profiles excluded.'},
        {'gse':'GSE107644','samples_in_GEO':106,'admitted_RNA_profiles':0,'assessable_broad_labels':0,'decision':'Inspected public barcode information has raw_file/processed_file/cell/barcode only; no explicit fine stage labels in that file. Atlas expression deferred; no stage inferred from cell names.'},
        {'gse':'GSE109033','samples_in_GEO':8,'admitted_RNA_profiles':0,'assessable_broad_labels':0,'decision':'GEO sample labels identify enrichments, not per-cell leptotene. Annotated Loupe data mentioned by paper require separate audit; no atlas ingested.'}]
    pd.DataFrame(review).to_csv(out/'candidate_review.csv',index=False)
    fixed={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in refpaths+[str(mapping_path.relative_to(ROOT)),str(paper.relative_to(ROOT))]}
    design={'scope':'Two new study sources; 6 stage-separated and 4 mixed-LZ/pachytene RNA profiles, 8 broad-label assessable and 2 preleptotene unscored.',
        'source_paper':provenance,'source_metadata_counts':{g:len(metadata(g)[1]) for g in ['GSE155967','GSE107644','GSE109033','GSE83264']},
        'GSE107644_barcode_rows':len(barcode),'GSE107644_barcode_columns':barcode.columns.tolist(),
        'genes':len(genes),'gene_folds_seed':20261024,'identifier_conversion':'Exact unversioned mouse Ensembl IDs to unique mouse symbols from previously frozen author genes.tsv; no online ID conversion, no effect-based alias choice.',
        'mapping_scope':'Old annotation may omit newer IDs; all exclusions and expression-mass coverage recorded. Direct source symbols in GSE83264 retained without alias remapping.',
        'normalization':'Same common-gene columns normalized before five gene folds; held-out genes contribute to totals. No held-out-specific scaling.',
        'reference_training':'No new reference trained or markers retuned. Three old reference hashes checked against the pre-existing frozen design.',
        'technical_caveat':'Full-length bulk HTSeq counts have gene-length/library-preparation effects unlike 3-prime UMI references; no lengths available for a validated correction. Compatibility results do not isolate biological stage mismatch.',
        'interpretation':'RNA weights are not cell-number fractions; broad label matches are not calibrated accuracy or proof of 100% purity. Mixed LZ is not a pure leptotene truth set.',
        'replication':'Two source-reported biological samples per stage; no inferential animal p-values, no gene-as-animal tests, no inference of exact animal counts from sample count.',
        'fixed_inputs_sha256':fixed,'not_claimed':['Genotype-verified WT status for every GSE155967 RNA sample','Calibrated cell proportions','A newly validated leptotene reference','Systematic literature coverage','Scite verification after plugin transport failure']}
    (out/'design.json').write_text(json.dumps(design,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    labels=[f'{r.gse}\n{r.stage}/{r.profile}' for r in evidence.itertuples()]
    for j,model in enumerate(MODELS):
        table=full[full.model==model].set_index('profile').loc[evidence.profile];m=metrics[metrics.model==model].set_index('profile').loc[evidence.profile]
        axes[0].scatter(np.arange(len(labels))+(j-1)*.15,table.Meiotic,s=30,label=model)
        axes[1].scatter(np.arange(len(labels))+(j-1)*.15,m.holdout_log2_RMSE,s=30,label=model)
    for ax in axes:ax.set(xticks=range(len(labels)),xticklabels=labels);ax.tick_params(axis='x',rotation=70,labelsize=7);ax.legend(fontsize=8)
    axes[0].set(ylabel='Fitted meiotic RNA weight (not cell fraction)',ylim=(0,1));axes[1].set(ylabel='Held-out gene log2 RMSE')
    fig.suptitle('Independent source compatibility; preleptotene displayed but not broad-label scored')
    fig.savefig(out/'stage_diagnostics.png',dpi=170);plt.close(fig)
    print(full[['gse','stage','profile','model','largest_class','Meiotic','matches_author_broad']].to_string(index=False),flush=True)


if __name__=='__main__':main()
