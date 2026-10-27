#!/usr/bin/env python3
"""Reconcile original Chen UMI cells with author stage and WT-QC workbooks."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from fetch_chen_umi import author_labels

STAGES=['A1','In','BS','BG2','G1','ePL','mPL','lPL','L','Z','eP','mP','lP','D','MI','MII','RS2','RS4','RS6','RS8']


def broad_stage(stage):
    if stage in ['A1','In','BS','BG2']:return 'Spermatogonia'
    if stage in ['G1','ePL','mPL','lPL']:return None
    if stage in ['L','Z','eP','mP','lP','D','MI','MII']:return 'Meiotic'
    if stage in ['RS2','RS4','RS6','RS8']:return 'Postmeiotic'
    raise ValueError('Unknown author stage')


def read_umi(path):
    x=pd.read_csv(path,sep='\t',index_col=0)
    a=x.to_numpy(dtype=float)
    if x.index.has_duplicates or x.columns.has_duplicates or x.index.isna().any() or not np.isfinite(a).all() or (a<0).any() or not (a==np.floor(a)).all():
        raise ValueError('Invalid source integer UMI matrix')
    return x


def main():
    out=ROOT/'results/reference_limits';out.mkdir(parents=True,exist_ok=True)
    spec=json.loads((ROOT/'data/chen_umi_assets.json').read_text())
    if spec['status']!='complete':raise ValueError('Source retrieval incomplete')
    for r in spec['assets']:
        path=ROOT/r['path']
        if path.is_symlink() or path.stat().st_size!=r['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest()!=r['sha256']:raise ValueError('Chen source cache changed')
    labels=author_labels(ROOT).set_index('matrix_cell')
    if set(labels.type)!=set(STAGES):raise ValueError('Author stage inventory differs')
    windows=[r for r in spec['assets'] if 'library' in r];common=None;seen=set();qc=[];inventory=[]
    assignments={}
    for r in windows:
        x=read_umi(ROOT/r['path'])
        if seen.intersection(x.columns):raise ValueError('Duplicate matrix cell across libraries')
        seen.update(x.columns);common=x.index if common is None else common.intersection(x.index)
        included=x.columns.intersection(labels.index)
        # Preserve the reviewed literal I/l cell-ID alias and independently
        # verify its QC fingerprint; stage always comes from author Table S2.
        if len(included):
            expected=labels.loc[included]
            totals=x[included].sum();detected=x[included].gt(0).sum()
            if not np.array_equal(totals.to_numpy(),expected['UMI'].to_numpy()) or not np.array_equal(detected.to_numpy(),expected['Gene Number'].to_numpy()):
                raise ValueError('Author WT QC UMI/detected-gene values do not reconcile')
            for cell in included:
                row=labels.loc[cell]
                qc.append({'cell':cell,'author_cell_id':row.cell,'explicit_ID_alias_used':bool(row.explicit_ID_alias_used),
                    'library':r['library'],'gsm':r['gsm'],'author_type':row['type'],'author_cluster':row['cluster'],
                    'original_barcode_join_present':bool(row.original_barcode_join_present),'verified_total_UMI':int(totals.loc[cell]),
                    'verified_detected_genes':int(detected.loc[cell]),'identity_method':'Matrix cell ID + author S2 + WT QC; explicit In/ln library+suffix alias only, exact totals/detected required; no expression-based stage inference.'})
        assignments[r['path']]=included.tolist()
        inventory.append({'library':r['library'],'gsm':r['gsm'],'matrix_cells':x.shape[1],'matrix_genes':x.shape[0],
            'author_stage_cells':len(included),'matrix_cells_without_author_S2_label':x.shape[1]-len(included)})
        print('Reconciled',r['library'],len(included),'labelled cells',flush=True)
    qc=pd.DataFrame(qc);missing=labels.loc[~labels.index.isin(qc.cell)].copy()
    # Never hide unresolved identifiers by retaining the expected stage counts.
    missing.to_csv(out/'chen_unresolved_cells.csv',index_label='cell')
    if len(missing):
        qc.to_csv(out/'chen_cell_QC.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        raise ValueError(f'{len(missing)} author stage cells absent by exact matrix ID; audit retained, no candidate reference generated')
    common=common.sort_values();totals={stage:np.zeros(len(common)) for stage in STAGES};means={stage:np.zeros(len(common)) for stage in STAGES}
    batch_sums={};batch_means={};batch_n={}
    counts=qc.author_type.value_counts()
    for r in windows:
        x=read_umi(ROOT/r['path'])
        for cell in assignments[r['path']]:
            stage=labels.loc[cell,'type'];v=x.loc[common,cell].to_numpy(dtype=float)
            totals[stage]+=v;means[stage]+=v/float(labels.loc[cell,'UMI'])*1e6
            if stage in ['L','Z','BS']:
                key=(stage,r['library']);batch_sums.setdefault(key,np.zeros(len(common)));batch_means.setdefault(key,np.zeros(len(common)))
                batch_sums[key]+=v;batch_means[key]+=v/float(labels.loc[cell,'UMI'])*1e6;batch_n[key]=batch_n.get(key,0)+1
    profiles=[];evidence=[]
    for weighting,table in [('RNA_sum',totals),('equal_cell_mean',means)]:
        for stage in STAGES:
            v=table[stage] if weighting=='RNA_sum' else table[stage]/counts[stage]
            name=f'Chen:{stage}:{weighting}';profiles.append(pd.Series(v,index=common,name=name))
            selected=qc[qc.author_type==stage]
            evidence.append({'profile':name,'author_type':stage,'weighting':weighting,'cells':int(counts[stage]),
                'source_libraries':int(selected.library.nunique()),'expected_broad_class':broad_stage(stage),
                'preparation':'WIN/RA synchronized fluorescent-marker germ cells, WT-labelled sheet; captured cells not in-vivo fractions',
                'author_stage_source':'Original paper Table S2; no cell-name or expression-based stage inference'})
    result=pd.concat(profiles,axis=1);result.index.name='gene';result.to_csv(out/'chen_stage_profiles.csv.gz',compression={'method':'gzip','mtime':0})
    batch_profiles=[];batch_evidence=[]
    for weighting,table in [('RNA_sum',batch_sums),('equal_cell_mean',batch_means)]:
        for (stage,library),values in sorted(table.items()):
            name=f'Chen:{stage}:{library}:{weighting}'
            v=values if weighting=='RNA_sum' else values/batch_n[(stage,library)]
            batch_profiles.append(pd.Series(v,index=common,name=name));batch_evidence.append({'profile':name,'author_type':stage,
                'library':library,'weighting':weighting,'cells':batch_n[(stage,library)]})
    batches=pd.concat(batch_profiles,axis=1);batches.index.name='gene';batches.to_csv(out/'chen_batch_prototypes.csv.gz',compression={'method':'gzip','mtime':0})
    pd.DataFrame(batch_evidence).to_csv(out/'chen_batch_evidence.csv',index=False)
    for stage in ['L','Z','BS']:
        keys=[k for k in batch_n if k[0]==stage]
        np.testing.assert_allclose(sum(batch_sums[k] for k in keys),result[f'Chen:{stage}:RNA_sum'])
        np.testing.assert_allclose(sum(batch_means[k] for k in keys)/sum(batch_n[k] for k in keys),result[f'Chen:{stage}:equal_cell_mean'])
    pd.DataFrame(evidence).to_csv(out/'chen_profile_evidence.csv',index=False)
    qc.to_csv(out/'chen_cell_QC.csv.gz',index=False,compression={'method':'gzip','mtime':0});pd.DataFrame(inventory).to_csv(out/'chen_library_inventory.csv',index=False)
    qc.groupby('author_type').agg(cells=('cell','size'),libraries=('library','nunique'),
        median_UMI=('verified_total_UMI','median'),median_detected_genes=('verified_detected_genes','median')).reindex(STAGES).to_csv(out/'chen_stage_inventory.csv')
    audit={'scope':'Complete author S2 stage/QC-selected set, not all experiments or animals in GSE107644.',
        'author_selected_cells':len(labels),'verified_matrix_cells':len(qc),'all_selected_cell_UMI_and_detected_gene_QC_match':True,
        'source_libraries':len(windows),'source_matrix_cells':int(sum(r['matrix_cells'] for r in inventory)),
        'all_library_common_gene_symbols':len(common),'source_matrix_gene_counts':sorted(set(r['matrix_genes'] for r in inventory)),
        'author_stages':len(STAGES),'profiles':len(evidence),'explicit_source_cell_ID_aliases_verified_with_QC':int(qc.explicit_ID_alias_used.sum()),
        'mutation_exclusion':'Only author S2 cells intersected with explicit Wild-type QC sheet; Sox30-cKO and Spo11-cKO sheets are not used.',
        'mapping':'Source unique symbols retained without alias guessing; all library-common genes used, no missing-feature zero imputation.',
        'limitations':['Stage-sorted synchronized transgenic-marker cells, not unbiased normal-tissue sampling.',
            'Source libraries and captured cells are not verified independent animals. Preleptotene is not broad-label scored.',
            'Reference construction and transfer below are exploratory follow-ups informed by prior mismatch, not untouched validation.']}
    (out/'chen_preparation_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':main()
