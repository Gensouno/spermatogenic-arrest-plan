#!/usr/bin/env python3
"""Stream author-labelled count triplets into library/type profiles, never dense cells."""
import json
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SOMATIC={'Endothelial_cells','Fetal_Leydig_1','Fetal_Leydig_2','Interstitial_tMg','Leydig_1','Leydig_2','PTM','Sertoli'}
MEIOTIC={'eP1','eP2','lP1','lP2','mP','D','MI','MII'}
WT_SAMPLES={'B6','P5','P10','P15','P20','P25','P30','P35'}


def broad_label(label):
    if label in SOMATIC:return 'Somatic'
    if label in MEIOTIC:return 'Meiotic'
    if label=='Spermatogonia':return 'Spermatogonia'
    if re.fullmatch(r'S(?:[1-9]|10|11)',label):return 'Postmeiotic'
    if label=='Outliers':return None
    raise ValueError('Unknown author annotation: '+label)


def parse_triplets(block, genes, cells):
    if not re.fullmatch(rb'[0-9 \t\r\n]+',block):raise ValueError('Invalid triplet characters')
    lines=block.splitlines()
    if any(len(line.split())!=3 for line in lines):raise ValueError('Expected exactly three fields per triplet row')
    values=np.fromstring(block.decode('ascii'),sep=' ',dtype=np.int64)
    if len(values)%3 or len(values)//3!=len(lines):raise ValueError('Malformed triplet rows')
    a=values.reshape(-1,3)
    if (a<1).any() or (a[:,0]>genes).any() or (a[:,1]>cells).any():raise ValueError('Invalid triplet value/index')
    key=a[:,1]*genes+a[:,0]
    if (np.diff(key)<=0).any():raise ValueError('Triplets must be unique and column-major sorted')
    return a


def source_metadata(base):
    cells=pd.read_csv(base/'cell_metadata.txt',sep=r'\s+')
    if cells.index.tolist()!=list(range(1,len(cells)+1)):raise ValueError('Author cell row order changed')
    if cells[['Library','Barcode']].duplicated().any():raise ValueError('Duplicate library/barcode')
    cells['expected_broad_class']=[broad_label(s) for s in cells.AnnotatedClusters]
    sdrf=pd.read_csv(base/'E-MTAB-6946.sdrf.txt',sep='\t')
    cols=['Source Name','Characteristics[age]','Characteristics[strain]','Characteristics[genotype]','Comment[cell count]']
    source=sdrf[cols].drop_duplicates().rename(columns=dict(zip(cols,['Library','age_days','strain','genotype','declared_cells'])))
    if source.Library.duplicated().any() or set(source.Library)!=set(cells.Library):raise ValueError('Ambiguous library/source mapping')
    observed=cells.groupby('Library').size()
    if not all(observed[r.Library]==r.declared_cells for r in source.itertuples()):raise ValueError('SDRF cell totals differ')
    source=source.merge(cells.groupby('Library').Sample.agg(['first','nunique']),left_on='Library',right_index=True)
    if (source['nunique']!=1).any():raise ValueError('Multiple sample names in library')
    source=source.rename(columns={'first':'Sample'}).drop(columns='nunique')
    source['normal_B6']=(source.genotype=='wild type genotype')&(source.strain=='C57BL6J')
    if set(source.loc[source.normal_B6,'Sample'])!=WT_SAMPLES:raise ValueError('Unexpected WT inventory')
    cells['normal_B6']=cells.Library.isin(source.loc[source.normal_B6,'Library'])
    return cells,source


def main():
    started=time.monotonic();base=ROOT/'data/raw/labeled_reference';out=ROOT/'results/labeled_reference'
    out.mkdir(parents=True,exist_ok=True)
    assets=json.loads((ROOT/'data/labeled_reference_assets.json').read_text())
    windows=[r for r in assets['assets'] if 'range' in r]
    selected_libraries={r['range']['library'] for r in windows}
    if selected_libraries!={'do17815','do26386','do18195'}:raise ValueError('Prespecified library subset differs')
    cells,source=source_metadata(base)
    if not source.set_index('Library').loc[sorted(selected_libraries),'normal_B6'].all():
        raise ValueError('A selected source library is not normal C57BL6J wild type')
    genes=pd.read_csv(base/'genes.tsv',sep='\t')
    if genes.ID.duplicated().any() or genes.ID.isna().any():raise ValueError('Duplicate/missing gene IDs')
    mouse=genes.ID.str.fullmatch(r'ENSMUSG[0-9]+')
    if not (mouse | genes.ID.str.fullmatch(r'ENSG[0-9]+')).all():raise ValueError('Unexpected feature species/ID')
    unique=mouse & genes.Symbol.notna() & ~genes.loc[mouse,'Symbol'].duplicated(False).reindex(genes.index,fill_value=False)
    use=cells.normal_B6 & cells.expected_broad_class.notna()
    groups=cells[use].groupby(['Sample','Library','AnnotatedClusters','expected_broad_class']).size().reset_index(name='cells')
    groups['selected_library']=groups.Library.isin(selected_libraries)
    groups['included']=(groups.cells>=20)&groups.selected_library
    groups=groups.merge(source[['Library','age_days']],on='Library',validate='many_to_one')
    groups.to_csv(out/'all_group_inventory.csv',index=False)
    chosen=groups[groups.included].reset_index(drop=True)
    lookup={(r.Library,r.AnnotatedClusters):i for i,r in chosen.iterrows()}
    assignment=np.array([lookup.get((r.Library,r.AnnotatedClusters),-1) if r.normal_B6 else -1 for r in cells.itertuples()])
    ng,nc,nk=len(genes),len(cells),len(chosen)
    sums=np.zeros((ng,nk));means=np.zeros_like(sums)
    totals=np.zeros(nc,dtype=np.int64);detected=np.zeros(nc,dtype=np.int64)
    denominators=cells.total_counts.to_numpy(dtype=np.float64)
    if not np.isfinite(denominators).all() or (denominators<=0).any():raise ValueError('Invalid author totals')
    processed=0;last=0
    rows,columns,full_nnz=map(int,assets['source_matrix']['header'][1].split())
    if (rows,columns)!=(ng,nc):raise ValueError('Metadata/matrix dimensions differ')
    nonzeros=sum(w['range']['nonzeros'] for w in windows)
    for window in windows:
        start_count=processed
        with (ROOT/window['path']).open('rb') as stream:
            while True:
                block=stream.read(4*1024**2)
                if not block:break
                if not block.endswith(b'\n'):block+=stream.readline()
                a=parse_triplets(block,ng,nc)
                if a[:,1].min()<window['range']['first_column'] or a[:,1].max()>window['range']['last_column']:
                    raise ValueError('Column outside declared library window')
                if a[0,1]*ng+a[0,0]<=last:raise ValueError('Unsorted/duplicate across chunks')
                last=int(a[-1,1]*ng+a[-1,0]);processed+=len(a)
                gi=a[:,0]-1;ci=a[:,1]-1;value=a[:,2]
                np.add.at(totals,ci,value);np.add.at(detected,ci,1)
                group=assignment[ci];keep=group>=0
                np.add.at(sums,(gi[keep],group[keep]),value[keep])
                np.add.at(means,(gi[keep],group[keep]),value[keep]/denominators[ci[keep]]*1e6)
        if processed-start_count!=window['range']['nonzeros']:raise ValueError('Window nonzero count differs')
        print('Aggregated',window['range']['library'],processed,'of',nonzeros,'entries;',round(time.monotonic()-started,1),'seconds',flush=True)
    if processed!=nonzeros:raise ValueError('Nonzero count differs from header')
    selected=cells.Library.isin(selected_libraries).to_numpy()
    if not np.array_equal(totals[selected],denominators[selected]):raise ValueError('Count sums do not match author per-cell totals')
    if not np.array_equal(detected[selected],cells.total_features_by_counts.to_numpy()[selected]):raise ValueError('Detected gene counts do not match author metadata')
    qc=cells.loc[selected,['Sample','Barcode','Library','AnnotatedClusters','normal_B6']].copy()
    qc['verified_count_total']=totals[selected];qc['verified_detected_genes']=detected[selected];qc['included_profile']=assignment[selected]>=0
    qc.to_csv(out/'cell_QC.csv.gz',index=False)
    profiles=[];evidence=[]
    for weighting,matrix in [('RNA_sum',sums),('equal_cell_mean',means/chosen.cells.to_numpy())]:
        x=matrix[unique]
        for i,row in chosen.iterrows():
            name=f'{row.Library}:{row.AnnotatedClusters}:{weighting}'
            evidence.append({'profile':name,**row.to_dict(),'weighting':weighting,
                'unique_mouse_expression_mass':float(x[:,i].sum()/matrix[:,i].sum())})
        profiles.append(pd.DataFrame(x,index=genes.loc[unique,'Symbol'],columns=[e['profile'] for e in evidence[-nk:]]))
    result=pd.concat(profiles,axis=1);result.index.name='gene'
    result.to_csv(out/'type_profiles.csv.gz')
    pd.DataFrame(evidence).to_csv(out/'profile_evidence.csv',index=False)
    source.to_csv(out/'library_evidence.csv',index=False)
    audit={'total_cells':nc,'normal_B6_cells':int(cells.normal_B6.sum()),'excluded_Tc_cells':int((~cells.normal_B6).sum()),
        'matrix_genes':ng,'full_matrix_nonzeros_declared':full_nnz,'selected_window_nonzeros':processed,
        'selected_libraries':sorted(selected_libraries),'selected_cells':int(selected.sum()),
        'all_selected_cell_totals_and_detected_genes_match':True,
        'mouse_features':int(mouse.sum()),'human_features_excluded':int((~mouse).sum()),'unique_mouse_symbols':int(unique.sum()),
        'normal_libraries':int(source.normal_B6.sum()),'normal_labels_excluding_outliers':int(groups.AnnotatedClusters.nunique()),
        'all_normal_nonoutlier_groups':len(groups),'included_groups_min20':nk,'profiles':result.shape[1],
        'included_cells':int((assignment>=0).sum()),'author_outliers_excluded':int((selected & cells.AnnotatedClusters.eq('Outliers')).sum()),
        'small_group_cells_excluded':int(groups.loc[groups.selected_library&~groups.included,'cells'].sum()),
        'label_provenance':'Original AnnotatedClusters supplied by authors; broad crosswalk fixed before expression fitting.',
        'order_evidence':'Author README binds metadata rows to matrix columns and gene table rows to matrix rows; all selected per-cell total and detected-feature counts independently reconciled. No full-matrix content hash claimed.',
        'elapsed_seconds':round(time.monotonic()-started,2)}
    (out/'preparation_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':main()
