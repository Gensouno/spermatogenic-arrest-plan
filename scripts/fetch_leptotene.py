#!/usr/bin/env python3
"""Retrieve six stage-separated RNA profiles and bounded comparator/label files."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import time

from discover_leptotene import ROOT,fetch

STAGES={'GSM5252960':'Preleptotene','GSM5252961':'Preleptotene','GSM5252962':'Leptotene',
    'GSM5252963':'Leptotene','GSM5252964':'Zygotene','GSM5252965':'Zygotene'}


def verify_asset(row,download):
    path=ROOT/row['path']
    if path.is_symlink() or not path.resolve().is_relative_to(ROOT.resolve()):raise ValueError('Unsafe source cache path')
    if not path.exists():
        if not download:raise FileNotFoundError('Missing frozen source '+row['path'])
        partial=path.with_name(path.name+'.restore-part')
        fetch(partial,row['url'],row['bytes'],True)
        if partial.stat().st_size!=row['bytes'] or hashlib.sha256(partial.read_bytes()).hexdigest()!=row['sha256']:
            raise ValueError('Restored source differs; unverified partial retained')
        partial.replace(path)
    if path.stat().st_size!=row['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
        raise ValueError('Frozen source differs: '+row['path'])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--download',action='store_true');args=parser.parse_args()
    raw=ROOT/'data/raw/leptotene';spec=ROOT/'data/leptotene_assets.json'
    if spec.exists():
        data=json.loads(spec.read_text())
        for row in data['assets']+[json.loads((ROOT/'data/leptotene_paper_asset.json').read_text())]:
            verify_asset(row,args.download)
        print('Verified',len(data['assets'])+1,'frozen leptotene inputs including primary paper',flush=True);return
    meta=ROOT/'results/leptotene_followup'
    discovery=json.loads((meta/'discovery.json').read_text());assets=list(discovery['assets']);selections=[]
    if discovery['metadata_retrieved']!=4:raise ValueError('Candidate metadata incomplete')
    d=json.loads((meta/'GSE155967_metadata.json').read_text())
    actual={gsm for gsm,s in d['samples'].items() if s.get('Sample_library_strategy')==['RNA-Seq']}
    if actual!=set(STAGES):raise ValueError('RNA sample inventory differs')
    for gsm,stage in STAGES.items():
        sample=d['samples'][gsm]
        if sample.get('Sample_organism_ch1')!=['Mus musculus'] or 'cell type: '+stage not in sample['Sample_characteristics_ch1']:
            raise ValueError('Species/stage identity differs')
        urls=[u for key,value in sample.items() if key.startswith('Sample_supplementary_file') for u in value]
        if len(urls)!=1 or not urls[0].endswith('.counts.txt.gz'):raise ValueError('Ambiguous RNA count file')
        url=urls[0].replace('ftp://ftp.ncbi.nlm.nih.gov/','https://ftp.ncbi.nlm.nih.gov/',1)
        row=fetch(raw/'expression'/Path(url).name,url,2*1024**2,args.download);assets.append(row)
        selections.append({'accession':'GSE155967','gsm':gsm,'stage':stage,'expected_broad_class':'Meiotic' if stage!='Preleptotene' else None,
            'path':row['path'],'original_title':sample['Sample_title'][0],
            'characteristics':sample['Sample_characteristics_ch1'],'genotype_evidence':'No per-sample genotype field; primary paper source/preparation must be audited before interpretation.'})
        print('Retrieved',gsm,stage,row['bytes'],'bytes',flush=True)
        time.sleep(.4)
    for gse,name in [('GSE83264','GSE83264_mRNA_genes_merged_counts.txt.gz'),
        ('GSE107644','GSE107644_barcode_information.txt.gz')]:
        d=json.loads((meta/(gse+'_metadata.json')).read_text());urls=d['series']['Series_supplementary_file']
        url=next(u for u in urls if u.endswith('/'+name)).replace('ftp://ftp.ncbi.nlm.nih.gov/','https://ftp.ncbi.nlm.nih.gov/',1)
        assets.append(fetch(raw/'expression'/name,url,8*1024**2,args.download));print('Retrieved',name,flush=True)
        time.sleep(.4)
    report={'retrieved_utc':datetime.now(timezone.utc).isoformat(),'scope':'Four candidate metadata files, six stage-separated RNA count files, one mixed-LZ comparator matrix, one barcode annotation file. No single-cell atlas expression ingested.',
        'selections':selections,'assets':assets,'count_reconciliation':{'candidate_series':4,'GSE155967_RNA_seq_samples':6,'selected_RNA_count_files':6}}
    spec.write_text(json.dumps(report,indent=2)+'\n')
    print('Frozen',len(assets),'source assets',flush=True)


if __name__=='__main__':main()
