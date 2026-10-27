#!/usr/bin/env python3
"""Reconcile complete P15 counts, then aggregate fixed author annotations."""
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from prepare_labeled_reference import parse_triplets,source_metadata

ROOT=Path(__file__).resolve().parents[1]
SOMATIC={'Immature_Leydig','Leydig_1','Leydig_2','Sertoli_1','Sertoli_2','PTM','Interstitial_tMg','Endothelial_cells'}
SG={'Differentiating_Spermatogonia_1','Differentiating_Spermatogonia_2','Undifferentiated_Spermatogonia'}
MEIOTIC={'Leptotene','Zygotene','Pachytene_1','Pachytene_2'}


def broad_label(label):
    if label in SOMATIC:return 'Somatic'
    if label in SG:return 'Spermatogonia'
    if label in MEIOTIC:return 'Meiotic'
    if label=='Outliers':return None
    raise ValueError('Unrecognized P15 author label: '+str(label))


def unique_mouse(genes):
    if genes.ID.isna().any() or genes.ID.duplicated().any():raise ValueError('Missing/duplicate feature ID')
    mouse=genes.ID.str.fullmatch(r'ENSMUSG[0-9]+')
    if not (mouse|genes.ID.str.fullmatch(r'ENSG[0-9]+')).all():raise ValueError('Unexpected species/ID')
    return mouse & genes.Symbol.notna() & ~genes.loc[mouse,'Symbol'].duplicated(False).reindex(genes.index,fill_value=False)


def mark_overlap(cells,old):
    keys=pd.MultiIndex.from_frame(old[['Library','Barcode']])
    if keys.duplicated().any() or cells[['Library','Barcode']].duplicated().any():raise ValueError('Duplicate library/barcode')
    return pd.MultiIndex.from_frame(cells[['Library','Barcode']]).isin(keys)


def verify_cells(totals,detected,cells):
    if not np.array_equal(totals,cells.total_counts.to_numpy()):raise ValueError('Per-cell count total differs')
    if not np.array_equal(detected,cells.total_features_by_counts.to_numpy()):raise ValueError('Per-cell detected features differ')


def aggregate(path,genes,cells,first,last,expected_nnz,groups,assignments):
    """One streaming pass; all cells checked even when excluded from profiles."""
    ng=len(genes);nc=len(cells);nk=len(groups)
    sums=np.zeros((ng,nk));means=np.zeros_like(sums)
    totals=np.zeros(last-first+1,dtype=np.int64);detected=np.zeros_like(totals)
    denom=cells.total_counts.to_numpy(dtype=float)
    if not np.isfinite(denom).all() or (denom<=0).any():raise ValueError('Invalid author totals')
    processed=0;previous=0
    with path.open('rb') as stream:
        while True:
            block=stream.read(4*1024**2)
            if not block:break
            if not block.endswith(b'\n'):block+=stream.readline()
            a=parse_triplets(block,ng,nc)
            key=a[:,1]*ng+a[:,0]
            if key[0]<=previous:raise ValueError('Duplicate/unsorted across chunks')
            previous=int(key[-1])
            gi=a[:,0]-1;ci=a[:,1]-1;v=a[:,2];local=a[:,1]-first
            if local.min()<0 or local.max()>=len(totals):raise ValueError('Triplet outside selected library')
            np.add.at(totals,local,v);np.add.at(detected,local,1)
            for assignment in assignments:
                group=assignment[ci];keep=group>=0
                np.add.at(sums,(gi[keep],group[keep]),v[keep])
                np.add.at(means,(gi[keep],group[keep]),v[keep]/denom[ci[keep]]*1e6)
            processed+=len(a)
    if processed!=expected_nnz:raise ValueError('Library nonzero total differs')
    verify_cells(totals,detected,cells.iloc[first-1:last])
    return sums,means/groups.cells.to_numpy(),totals,detected


def main():
    begun=time.monotonic();out=ROOT/'results/emptydrops_p15';out.mkdir(parents=True,exist_ok=True)
    raw=ROOT/'data/raw/labeled_reference'
    cells=pd.read_csv(raw/'cell_metadata_emptyDrops.txt',sep=r'\s+',dtype={'P15Clusters':'string'},low_memory=False)
    if cells.index.tolist()!=list(range(1,len(cells)+1)):raise ValueError('Metadata order differs')
    old,source=source_metadata(raw)
    info=source.set_index('Library').loc['do18195']
    if not info.normal_B6 or info.Sample!='P15' or int(info.age_days)!=15:raise ValueError('Source age/strain/genotype differs')
    selected=cells.Library.eq('do18195')
    if int(selected.sum())!=13168 or not cells.loc[selected,'Sample'].eq('P15').all():raise ValueError('Selected library differs')
    cells['overlap']=mark_overlap(cells,old)
    cells['expected_broad_class']=None
    cells.loc[selected,'expected_broad_class']=[broad_label(s) for s in cells.loc[selected,'P15Clusters']]
    inventory=[];assignments=[]
    for stratum,mask in [('full',selected),('overlap',selected&cells.overlap),('EmptyDrops_only',selected&~cells.overlap)]:
        counts=cells[mask].groupby(['P15Clusters'],dropna=False).size()
        for label,n in counts.items():
            expected=broad_label(label)
            inventory.append({'branch':'EmptyDrops','stratum':stratum,'author_type':label,'cells':int(n),
                'expected_broad_class':expected,'included':expected is not None and n>=20,
                'exclusion_reason':'author_outlier' if expected is None else ('fewer_than_20_cells' if n<20 else '')})
    inv=pd.DataFrame(inventory);chosen=inv[inv.included].reset_index(drop=True)
    for stratum in ['full','overlap','EmptyDrops_only']:
        lookup={r.author_type:i for i,r in chosen[chosen.stratum==stratum].iterrows()}
        eligible=selected & (True if stratum=='full' else (cells.overlap if stratum=='overlap' else ~cells.overlap))
        assignments.append(np.array([lookup.get(label,-1) if ok else -1 for label,ok in zip(cells.P15Clusters,eligible)]))
    observed=sum((np.bincount(a[a>=0],minlength=len(chosen)) for a in assignments))
    if not np.array_equal(observed,chosen.cells.to_numpy()):raise ValueError('Profile assignment cell counts differ')
    genes=pd.read_csv(ROOT/'data/raw/emptydrops_p15/genes_emptyDrops.tsv',sep='\t')
    unique=unique_mouse(genes)
    assets=json.loads((ROOT/'data/emptydrops_p15_assets.json').read_text())
    window=next(r for r in assets['assets'] if 'range' in r)
    w=window['range']
    indices=cells.index[selected].to_numpy()
    if (not np.array_equal(indices,np.arange(indices[0],indices[-1]+1)) or
            (w['first_column'],w['last_column'],w['cells'],w['nonzeros'])!=
            (int(indices[0]),int(indices[-1]),len(indices),int(cells.loc[selected,'total_features_by_counts'].sum()))):
        raise ValueError('Frozen range does not describe the complete selected library')
    if tuple(map(int,assets['source_matrix']['header'][1].split()))!=(len(genes),len(cells),int(cells.total_features_by_counts.sum())):
        raise ValueError('Source matrix dimensions differ from all gene/cell metadata')
    sums,means,totals,detected=aggregate(ROOT/window['path'],genes,cells,w['first_column'],w['last_column'],w['nonzeros'],chosen,assignments)
    print('Verified and aggregated all 13,168 EmptyDrops cells',flush=True)
    qc=cells.loc[selected,['Library','Barcode','P15Clusters','overlap','expected_broad_class']].copy()
    qc['verified_count_total']=totals;qc['verified_detected_genes']=detected
    qc.to_csv(out/'cell_QC.csv.gz',index=False)
    # Reaggregate the SAME barcodes from the CellRanger count branch using the
    # same EmptyDrops P15 labels. This is a controlled count-branch comparison.
    oldgenes=pd.read_csv(raw/'genes.tsv',sep='\t');oldunique=unique_mouse(oldgenes)
    shared=cells.loc[selected&cells.overlap,['Library','Barcode','P15Clusters']]
    labels=old[['Library','Barcode']].merge(shared,on=['Library','Barcode'],how='left',validate='one_to_one').P15Clusters
    oldchosen=chosen[chosen.stratum=='overlap'].copy().reset_index(drop=True);oldchosen['branch']='CellRanger'
    lookup={r.author_type:i for i,r in oldchosen.iterrows()}
    assignment=np.array([lookup.get(s,-1) for s in labels])
    if not np.array_equal(np.bincount(assignment[assignment>=0],minlength=len(oldchosen)),oldchosen.cells.to_numpy()):
        raise ValueError('Matched barcode profile cell counts differ')
    oldassets=json.loads((ROOT/'data/labeled_reference_assets.json').read_text())
    oldwindow=next(r for r in oldassets['assets'] if r.get('range',{}).get('library')=='do18195');ow=oldwindow['range']
    oldsum,oldmean,oldtotal,olddetected=aggregate(ROOT/oldwindow['path'],oldgenes,old,ow['first_column'],ow['last_column'],ow['nonzeros'],oldchosen,[assignment])
    oldqc=old.loc[old.Library.eq('do18195'),['Library','Barcode','AnnotatedClusters']].copy()
    oldqc['P15Clusters']=labels.to_numpy()[old.Library.eq('do18195').to_numpy()]
    oldqc['verified_count_total']=oldtotal;oldqc['verified_detected_genes']=olddetected
    oldqc.to_csv(out/'cellranger_cell_QC.csv.gz',index=False)
    def frames(gene_table,keep,group_table,matrices):
        frames=[];rows=[]
        for weighting,matrix in matrices:
            x=matrix[keep]
            for i,row in group_table.iterrows():
                name=f'{row.branch}:{row.stratum}:{row.author_type}:{weighting}'
                rows.append({'profile':name,**row.to_dict(),'weighting':weighting,
                    'unique_mouse_expression_mass':float(x[:,i].sum()/matrix[:,i].sum())})
            frames.append(pd.DataFrame(x,index=gene_table.loc[keep,'Symbol'],columns=[r['profile'] for r in rows[-len(group_table):]]))
        return pd.concat(frames,axis=1),rows
    x,ev=frames(genes,unique,chosen,[('RNA_sum',sums),('equal_cell_mean',means)])
    ox,oev=frames(oldgenes,oldunique,oldchosen,[('RNA_sum',oldsum),('equal_cell_mean',oldmean)])
    # The diagnostic uses this identical symbol universe for every new profile.
    common=x.index.intersection(ox.index).sort_values()
    count_comparison=[]
    for label in oldchosen.author_type:
        a=x.loc[common,f'EmptyDrops:overlap:{label}:RNA_sum'].to_numpy()
        b=ox.loc[common,f'CellRanger:overlap:{label}:RNA_sum'].to_numpy()
        count_comparison.append({'author_type':label,'genes':len(common),
            'different_gene_sums':int(np.count_nonzero(a!=b)),'max_absolute_gene_sum_difference':float(np.max(np.abs(a-b))),
            'EmptyDrops_common_gene_count_sum':float(a.sum()),'CellRanger_common_gene_count_sum':float(b.sum())})
    evidence=pd.DataFrame(ev+oev)
    combined=pd.concat([x.loc[common],ox.loc[common]],axis=1);combined.index.name='gene'
    for table in [x,ox]:
        for name in table:
            evidence.loc[evidence.profile.eq(name),'branch_common_expression_mass']=table.loc[common,name].sum()/table[name].sum()
    combined.to_csv(out/'type_profiles.csv.gz');evidence.to_csv(out/'profile_evidence.csv',index=False)
    pd.concat([inv,oldchosen],ignore_index=True).to_csv(out/'group_inventory.csv',index=False)
    qc.groupby(['P15Clusters','overlap']).agg(cells=('Barcode','size'),median_counts=('verified_count_total','median'),
        median_detected=('verified_detected_genes','median')).reset_index().to_csv(out/'cell_complexity_summary.csv',index=False)
    audit={'scope':'Complete do18195 P15 EmptyDrops library and matched CellRanger P15 count-branch comparison.',
        'selected_cells':int(selected.sum()),'overlap_cells':int((selected&cells.overlap).sum()),
        'emptydrops_only_cells':int((selected&~cells.overlap).sum()),'cellranger_P15_cells':int(old.Library.eq('do18195').sum()),
        'emptydrops_nonzeros':w['nonzeros'],'cellranger_nonzeros':ow['nonzeros'],
        'all_selected_cell_totals_and_detected_features_match':True,
        'emptydrops_genes':len(genes),'emptydrops_unique_mouse_symbols':int(unique.sum()),
        'cellranger_unique_mouse_symbols':int(oldunique.sum()),'branch_common_symbols':len(common),
        'matched_aggregate_count_comparison':count_comparison,
        'count_comparison_scope':'Gene-wise aggregated counts for identical eligible barcode groups; not a per-cell gene-count identity assertion.',
        'author_outliers_excluded':int(qc.P15Clusters.eq('Outliers').sum()),
        'included_profile_groups':len(chosen)+len(oldchosen),'weighting_views':2,'profiles':combined.shape[1],
        'source_age_days':int(info.age_days),'source_genotype':info.genotype,'source_strain':info.strain,
        'SDRF_cell_count_scope':'The SDRF count describes CellRanger, not EmptyDrops; not used as an EmptyDrops count assertion.',
        'label_provenance':'Unmodified author P15Clusters, matched by Library+Barcode. Same labels applied to both count branches.',
        'replication':'One library; strata overlap and weighting views are not independent biological replicates.',
        'elapsed_seconds':round(time.monotonic()-begun,2)}
    (out/'preparation_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':main()
