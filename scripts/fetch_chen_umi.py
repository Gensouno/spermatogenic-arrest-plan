#!/usr/bin/env python3
"""Fetch only author-labelled WT-sheet UMI libraries; no mutant or full-length data."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import time
import zipfile

import pandas as pd

from discover_leptotene import ROOT,fetch
from discover_expansion import parse_soft
import gzip

LABEL_FILE='41422_2018_74_MOESM26_ESM.xlsx'
QC_FILE='41422_2018_74_MOESM25_ESM.xlsx'


def restore_asset(row,download):
    path=ROOT/row['path']
    if path.is_symlink() or not path.resolve().is_relative_to(ROOT.resolve()):raise ValueError('Unsafe Chen source path')
    if not path.exists():
        partial=path.with_name(path.name+'.restore-part')
        if partial.is_symlink():raise ValueError('Symlink partial refused')
        if 'source_bundle' in row:
            with zipfile.ZipFile(ROOT/row['source_bundle']) as z:data=z.read(row['member'])
            if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:raise ValueError('Derived workbook differs from frozen source')
            partial.parent.mkdir(parents=True,exist_ok=True);partial.write_bytes(data)
        else:
            if not download:raise FileNotFoundError('Missing Chen source; use --download: '+row['path'])
            fetch(partial,row['url'],row['bytes'],True)
        if partial.stat().st_size!=row['bytes'] or hashlib.sha256(partial.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Restored Chen source differs; unverified partial retained')
        partial.replace(path)
    if path.stat().st_size!=row['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Existing Chen cache differs')


def reconcile_id_aliases(joined,barcode):
    joined=joined.copy()
    for col in ['raw_file','processed_file','barcode']:joined[col]=joined[col].astype('string')
    joined['original_barcode_join_present']=joined.raw_file.notna()
    joined['matrix_cell']=joined.cell;joined['explicit_ID_alias_used']=False
    for i,row in joined[joined.raw_file.isna()].iterrows():
        match=re.fullmatch(r'(SpSC(?:24|34))_In_sc\.([0-9]+)',row.cell)
        if match:
            alias=f'{match[1]}_ln_sc.{match[2]}'
            candidate=barcode[barcode.cell.eq(alias)]
            if len(candidate)!=1 or candidate.iloc[0].raw_file!=match[1]:raise ValueError('Unresolved explicit cell-ID alias')
            joined.loc[i,['raw_file','processed_file','barcode']]=candidate.iloc[0][['raw_file','processed_file','barcode']].to_numpy()
            joined.loc[i,'matrix_cell']=alias;joined.loc[i,'explicit_ID_alias_used']=True
    if joined.matrix_cell.duplicated().any():raise ValueError('Cell-ID alias collision')
    return joined


def author_labels(root):
    base=root/'data/raw/reference_limits/chen_labels'
    labels=pd.read_excel(base/(LABEL_FILE+'.verified'),sheet_name='Sheet1',header=4)
    qc=pd.read_excel(base/(QC_FILE+'.verified'),sheet_name='Wild-type QC',header=4).rename(columns={'Cell':'cell'})
    barcode=pd.read_csv(root/'data/raw/leptotene/expression/GSE107644_barcode_information.txt.gz',sep='\t')
    if labels.columns.tolist()!=['cell','type','cluster'] or len(labels)!=1136 or labels.cell.duplicated().any():raise ValueError('Unexpected author stage labels')
    if len(qc)!=1204 or qc.cell.duplicated().any() or barcode.cell.duplicated().any():raise ValueError('Unexpected WT QC/barcode inventory')
    joined=labels.merge(qc,on='cell',how='left',validate='one_to_one').merge(barcode,on='cell',how='left',validate='one_to_one')
    if joined['UMI'].isna().any():raise ValueError('Author stage cell missing in WT QC sheet')
    # Author S1/S2 write capital-I "In"; source matrix/barcode tables write
    # lowercase-L "ln" for these two libraries. Resolve only this literal pair,
    # requiring unique library+cell suffix and later exact UMI/detected QC.
    joined=reconcile_id_aliases(joined,barcode)
    # Prefix used only as a library identifier, never as an inferred stage.
    joined['library_candidate']=joined.raw_file.fillna(joined.cell.str.extract(r'^(SpSC[0-9]+)_',expand=False))
    if joined.library_candidate.isna().any():raise ValueError('Cell library unresolved')
    return joined


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--download',action='store_true');args=parser.parse_args()
    spec=ROOT/'data/chen_umi_assets.json'
    if spec.exists():
        old=json.loads(spec.read_text())
        for r in old['assets']:
            restore_asset(r,args.download)
        if old['status']=='complete':print('Verified',len(old['assets']),'Chen source assets; no network',flush=True);return
    labels=author_labels(ROOT)
    meta=ROOT/'data/raw/leptotene/metadata/GSE107644_family.soft.gz'
    _,samples=parse_soft(gzip.decompress(meta.read_bytes()).decode())
    available={}
    for gsm,r in samples.items():
        for key,urls in r.items():
            if not key.startswith('Sample_supplementary_file'):continue
            for u in urls:
                match=re.fullmatch(r'GSM[0-9]+_(SpSC[0-9]+)_UMIcount.txt.gz',Path(u).name)
                if match:
                    if match[1] in available:raise ValueError('Ambiguous UMI source library')
                    if r.get('Sample_organism_ch1')!=['Mus musculus'] or r.get('Sample_library_strategy')!=['RNA-Seq']:raise ValueError('Unexpected source organism/strategy')
                    available[match[1]]={'gsm':gsm,'url':u.replace('ftp://ftp.ncbi.nlm.nih.gov/','https://ftp.ncbi.nlm.nih.gov/',1),'library':match[1]}
    selected=sorted(set(labels.library_candidate))
    if not set(selected).issubset(available) or len(selected)>40:raise ValueError('Required WT-stage libraries missing/exceed budget')
    oldassets=json.loads((ROOT/'data/leptotene_assets.json').read_text())['assets']
    assets=[r for r in oldassets if Path(r['path']).name in ['GSE107644_family.soft.gz','GSE107644_barcode_information.txt.gz']]
    bundle=ROOT/'data/raw/reference_limits/chen_labels/PMC6123400_supplementary.zip'
    assets.append({'path':str(bundle.relative_to(ROOT)),'url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6123400/supplementaryFiles',
        'bytes':bundle.stat().st_size,'sha256':hashlib.sha256(bundle.read_bytes()).hexdigest()})
    with zipfile.ZipFile(bundle) as z:
        for name in [LABEL_FILE,QC_FILE]:
            member=next(n for n in z.namelist() if Path(n).name==name);p=bundle.parent/(name+'.verified')
            data=z.read(member)
            if p.is_symlink() or p.read_bytes()!=data or not data.startswith(b'PK'):raise ValueError('Workbook does not match source bundle')
            assets.append({'path':str(p.relative_to(ROOT)),'source_bundle':str(bundle.relative_to(ROOT)),'member':member,
                'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
    report={'status':'partial','retrieved_utc':datetime.now(timezone.utc).isoformat(),
        'scope':'Whole official UMI libraries containing author S2 stage-labelled cells in WT QC sheet, not all GSE107644 assays.',
        'author_S2_cells':len(labels),'direct_barcode_join_cells':int(labels.raw_file.notna().sum()),
        'unresolved_barcode_table_cells':int(labels.raw_file.isna().sum()),'selected_libraries':selected,'assets':assets}
    for library in selected:
        source=available[library];path=ROOT/'data/raw/reference_limits/chen_umi'/Path(source['url']).name
        r=fetch(path,source['url'],16*1024**2,args.download)
        with gzip.open(path,'rb') as stream:
            if b'\t' not in stream.readline():raise ValueError('Not a tab-separated UMI matrix')
        assets.append({**r,**source})
        if sum(r['bytes'] for r in assets)>400*1024**2:raise ValueError('Download resource cap exceeded')
        spec.write_text(json.dumps(report,indent=2)+'\n')
        print('Verified',library,len([r for r in assets if 'library' in r]),'of',len(selected),r['bytes'],'bytes',flush=True)
        time.sleep(.4)
    report['status']='complete';spec.write_text(json.dumps(report,indent=2)+'\n')
    print('Frozen',len(assets),'source assets;',len(selected),'UMI libraries',flush=True)


if __name__=='__main__':main()
