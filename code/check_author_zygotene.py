#!/usr/bin/env python3
"""Use explicit author P15 fine labels only for barcodes already in cached counts."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from prepare_labeled_reference import parse_triplets
from diagnose_reference_transfer import MODELS,broad_weights,gene_folds,train_fit

ROOT=Path(__file__).resolve().parents[1]


def join_author_labels(original,fine):
    keys=['Library','Barcode']
    if original[keys].duplicated().any() or fine[keys].duplicated().any():raise ValueError('Ambiguous author barcode join')
    old=original.copy();old['original_matrix_column']=old.index
    return old.merge(fine[keys+['P15Clusters']],on=keys,how='left',validate='one_to_one')


def main():
    base=ROOT/'data/raw/labeled_reference';out=ROOT/'results/labeled_reference'
    original=pd.read_csv(base/'cell_metadata.txt',sep=r'\s+')
    fine=pd.read_csv(base/'cell_metadata_emptyDrops.txt',sep=r'\s+',dtype={'P15Clusters':'string'},low_memory=False)
    old=original[original.Library=='do18195'];new=fine[fine.Library=='do18195']
    joined=join_author_labels(old,new)
    labels=joined.P15Clusters.value_counts()
    selected=joined[joined.P15Clusters=='Zygotene']
    if len(selected)<20:raise ValueError('Insufficient matched author Zygotene cells')
    mask=np.zeros(len(original),dtype=bool);mask[selected.original_matrix_column.to_numpy()-1]=True
    denominator=original.total_counts.to_numpy(dtype=float)
    genes=pd.read_csv(base/'genes.tsv',sep='\t');ng=len(genes)
    counts=np.zeros((ng,2));totals=np.zeros(len(original),dtype=np.int64);detected=np.zeros_like(totals)
    with (base/'do18195.triplets').open('rb') as stream:
        while True:
            block=stream.read(4*1024**2)
            if not block:break
            if not block.endswith(b'\n'):block+=stream.readline()
            a=parse_triplets(block,ng,len(original));gi=a[:,0]-1;ci=a[:,1]-1;v=a[:,2]
            keep=mask[ci]
            counts[:,0]+=np.bincount(gi[keep],weights=v[keep],minlength=ng)
            counts[:,1]+=np.bincount(gi[keep],weights=v[keep]/denominator[ci[keep]]*1e6,minlength=ng)/len(selected)
            np.add.at(totals,ci[keep],v[keep]);np.add.at(detected,ci[keep],1)
    if not np.array_equal(totals[mask],denominator[mask]) or not np.array_equal(detected[mask],original.total_features_by_counts.to_numpy()[mask]):
        raise ValueError('Matched Zygotene count/QC mismatch')
    mouse=genes.ID.str.fullmatch(r'ENSMUSG[0-9]+')
    unique=mouse & genes.Symbol.notna() & ~genes.loc[mouse,'Symbol'].duplicated(False).reindex(genes.index,fill_value=False)
    x=pd.DataFrame(counts[unique],index=genes.loc[unique,'Symbol'],columns=['RNA_sum','equal_cell_mean'])
    refs={m:pd.read_csv(ROOT/f'results/stage_reference/pooled_{m}_reference.csv.gz',index_col=0) for m in MODELS}
    foldmap=gene_folds(sorted(refs['centroid_4'].index));shared=refs['centroid_4'].index.intersection(x.index).sort_values()
    y=x.loc[shared].div(x.loc[shared].sum(),axis=1)*1e6
    refs={m:r.loc[shared].div(r.loc[shared].sum(),axis=1)*1e6 for m,r in refs.items()}
    scale=1/np.sqrt(refs['centroid_4'].mean(axis=1).to_numpy()+1);folds=np.array([foldmap[g] for g in shared])
    metrics=[];weights=[]
    for model,ref in refs.items():
        s=ref.to_numpy()
        for weighting in x:
            target=y[weighting].to_numpy();prediction=np.empty_like(target)
            for fold in [-1,0,1,2,3,4]:
                coef=train_fit(s,target,scale,np.ones(len(shared),dtype=bool) if fold==-1 else folds!=fold)
                if fold!=-1:prediction[folds==fold]=s[folds==fold]@coef
                p=broad_weights(coef,ref.columns,model)
                weights.append({'model':model,'weighting':weighting,'fold':fold,'cells':len(selected),
                    'largest_class':p.idxmax(),'expected_broad_class':'Meiotic','expected_class_weight':p['Meiotic'],**p.to_dict()})
            metrics.append({'model':model,'weighting':weighting,'cells':len(selected),'genes':len(shared),
                'holdout_log2_RMSE':float(np.sqrt(np.mean((np.log2(target+.5)-np.log2(prediction+.5))**2)))})
    pd.DataFrame(metrics).to_csv(out/'zygotene_overlap_metrics.csv',index=False)
    pd.DataFrame(weights).to_csv(out/'zygotene_overlap_weights.csv',index=False)
    audit={'scope':'Secondary author-label reassignment sensitivity on already cached CellRanger P15 expression; no EmptyDrops expression matrix processed.',
        'matched_P15_barcodes':int(joined.P15Clusters.notna().sum()),'matched_author_fine_label_counts':{str(k):int(v) for k,v in labels.items()},
        'Zygotene_cells_with_cached_expression':len(selected),'all_author_EmptyDrops_Zygotene_cells':int((new.P15Clusters=='Zygotene').sum()),
        'Leptotene_cells_with_cached_expression':int(labels.get('Leptotene',0)),
        'Leptotene_not_fitted':'Below prespecified >=20 cells/group threshold; do not merge with Zygotene to conceal coverage gap.',
        'Zygotene_original_AnnotatedClusters':selected.AnnotatedClusters.value_counts().to_dict(),
        'matched_cell_counts_and_detected_features_verified':True,
        'limitations':['The matched Zygotene subset excludes most EmptyDrops-only low-complexity cells and is not representative of all 938 author-labelled Zygotene cells.',
            'Author fine labels come from a different call/filter/annotation branch; this is an explicit sensitivity join, not replacement of the original primary labels.',
            'One library, two weighting views, no animal replication or in-vivo cell-fraction accuracy.']}
    (out/'zygotene_overlap_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('Author Zygotene overlap:',len(selected),'cells; Leptotene available:',int(labels.get('Leptotene',0)),flush=True)


if __name__=='__main__':main()
