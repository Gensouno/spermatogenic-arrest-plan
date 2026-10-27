#!/usr/bin/env python3
"""Inspect author EmptyDrops metadata only; do not infer L/Z from a generic NewCell tag."""
import argparse
import json
from pathlib import Path

import pandas as pd
from fetch_labeled_reference import CACHE,ROOT,fetch,file_entries


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--download',action='store_true');args=p.parse_args()
    name='cell_metadata_emptyDrops.txt'
    row=fetch(name,'https://www.ebi.ac.uk/biostudies/files/E-MTAB-6946/'+name,15*1024**2,args.download)
    study=json.loads((CACHE/'study.json').read_text());declared={r['path']:r for r in file_entries(study)}
    if row['bytes']!=declared[name]['size']:raise ValueError('Declared metadata size differs')
    path=ROOT/'data/emptydrops_label_assets.json'
    manifest={'scope':'One author EmptyDrops annotation file only; no EmptyDrops expression matrix downloaded or analysed.',
        'assets':[row]}
    if path.exists() and json.loads(path.read_text())!=manifest:raise ValueError('Frozen EmptyDrops metadata differs')
    if not path.exists():path.write_text(json.dumps(manifest,indent=2)+'\n')
    x=pd.read_csv(CACHE/name,sep=r'\s+',dtype={'P15Clusters':'string'},low_memory=False)
    if x[['Library','Barcode']].duplicated().any():raise ValueError('Duplicate library/barcode')
    p15=x[x.Sample=='P15']
    fields=[c for c in x if any(t in c.lower() for t in ['cluster','type','p15'])]
    original=pd.read_csv(CACHE/'cell_metadata.txt',sep=r'\s+')
    original_keys=set(zip(original.Library,original.Barcode))
    overlap=sum((r.Library,r.Barcode) in original_keys for r in x.itertuples())
    p15_overlap=sum((r.Library,r.Barcode) in original_keys for r in p15.itertuples())
    p15_original=int((original.Sample=='P15').sum())
    lz=p15[p15.P15Clusters.isin(['Leptotene','Zygotene'])]
    audit={'scope':manifest['scope'],'rows':len(x),'columns':x.columns.tolist(),'P15_rows':len(p15),
        'library_counts':x.groupby(['Sample','Library']).size().reset_index(name='cells').to_dict('records'),
        'annotation_columns':fields,
        'P15_annotation_values':{c:{str(k):int(v) for k,v in p15[c].value_counts(dropna=False).items()} for c in fields},
        'CellRanger_library_barcode_overlap':overlap,
        'P15_CellRanger_rows':p15_original,'P15_CellRanger_overlap':p15_overlap,
        'P15_EmptyDrops_only_barcodes':len(p15)-p15_overlap,'P15_CellRanger_only_barcodes':p15_original-p15_overlap,
        'explicit_P15_LZ_author_labels_present':{'Leptotene','Zygotene'}.issubset(set(p15.P15Clusters.dropna())),
        'P15_LZ_author_label_counts':{str(k):int(v) for k,v in lz.P15Clusters.value_counts().items()},
        'P15_LZ_versus_generic_annotations':pd.crosstab(lz.P15Clusters,lz.AnnotatedClusters).to_dict(),
        'warning':'Use explicit P15Clusters Leptotene/Zygotene, never relabel all NewCell as L/Z. Processed public CellRanger and EmptyDrops barcode sets are not strict subsets; do not infer gains from simple total-row subtraction.',
        'expression_matrix_processed':False}
    out=ROOT/'results/labeled_reference';out.mkdir(parents=True,exist_ok=True)
    (out/'emptydrops_label_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':main()
