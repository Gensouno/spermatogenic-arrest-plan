#!/usr/bin/env python3
"""Stream a developmental atlas into library-level diagnostic profiles."""
import csv
import gzip
import json
from pathlib import Path
import re
import time
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PREFIXES = ['PND06','PND14','PND18pre','PND18','PND25','PND30','Ad-Thy1','Ad-depleteSp']
GSMS = dict(zip(PREFIXES, [f'GSM{3449736+i}' for i in range(8)]))


def read_author_table(path):
    # This author workbook uses ISO strict OOXML, not the transitional namespace
    # accepted by openpyxl. Read only the small fixed data sheet, never macros.
    with zipfile.ZipFile(path) as archive:
        strings = ET.fromstring(archive.read('xl/sharedStrings.xml'))
        shared = [''.join(x.itertext()) for x in strings]
        sheet = ET.fromstring(archive.read('xl/worksheets/sheet1.xml'))
    rows = []
    for row in sheet.iter():
        if row.tag.rsplit('}',1)[-1] != 'row': continue
        values = []
        for cell in row:
            v = next((x.text for x in cell if x.tag.rsplit('}',1)[-1]=='v'),None)
            if v is None: raise ValueError('Missing sample table value')
            values.append(shared[int(v)] if cell.attrib.get('t')=='s' else v)
        if len(values)!=2: raise ValueError('Unexpected author table width')
        rows.append(values)
    if len(rows)!=9: raise ValueError('Expected header and eight author library rows')
    return pd.DataFrame(rows[1:],columns=rows[0])


def parse_count_row(line, n):
    name, separator, text = line.rstrip('\r\n').partition(',')
    if not separator: raise ValueError('Missing count fields')
    gene = next(csv.reader([name]))[0]
    # Counting separators avoids np.fromstring's silent truncation on bad tails.
    if text.count(',')!=n-1: raise ValueError('Wrong number of count fields')
    x = np.fromstring(text,sep=',',dtype=np.float64)
    if (len(x)!=n or not np.isfinite(x).all() or (x<0).any()
            or (np.abs(x-np.round(x))>1e-6).any()):
        raise ValueError('Invalid nonnegative integer counts')
    return gene,x


def barcode_assignments(columns):
    if len(columns)!=len(set(columns)): raise ValueError('Duplicate barcodes')
    observed=[]
    for barcode in columns:
        prefix,separator,sequence=barcode.rpartition('_')
        if not separator or prefix not in PREFIXES or not re.fullmatch('[ACGT]+',sequence):
            raise ValueError('Unrecognized source library/barcode format')
        observed.append(PREFIXES.index(prefix))
    if set(observed)!=set(range(8)): raise ValueError('Missing a source library')
    return np.array(observed,dtype=int)


def main():
    started=time.monotonic()
    base=ROOT/'data/raw/age_reference'
    path=base/'GSE121904_testis_maturation_scRNAseq_raw_mapped_counts.csv.gz'
    out=ROOT/'results/age_reference';out.mkdir(parents=True,exist_ok=True)
    meta=json.loads((out/'metadata.json').read_text())['samples']
    author=read_author_table(base/'author/doc/20180530_sample_info.xlsx')
    author.to_csv(out/'author_sample_table.csv',index=False)
    sample_names=author.iloc[:,0].tolist()
    if set(sample_names)!=set(PREFIXES): raise ValueError('Author table/prefix inventory differs')
    with gzip.open(path,'rt') as stream:
        header=next(csv.reader([stream.readline()]))
        columns=header[1:];assignment=barcode_assignments(columns)
        n=len(columns);totals=np.zeros(n);detected=np.zeros(n,dtype=int);mt=np.zeros(n)
        genes=[];raw=[]
        for line in stream:
            gene,x=parse_count_row(line,n);genes.append(gene)
            totals+=x;detected+=x>0
            if gene.startswith('mt-'):mt+=x
            raw.append(np.bincount(assignment,weights=x,minlength=8))
            if len(genes)%5000==0:print('PASS1',len(genes),'genes',round(time.monotonic()-started,1),'seconds',flush=True)
    if len(genes)!=len(set(genes)):raise ValueError('Duplicate gene symbols')
    if (totals<=0).any():raise ValueError('Empty library/barcode')
    mitochondrial=mt/totals
    paper_qc=(detected>500)&(detected<8000)&(totals>2000)&(totals<125000)&(mitochondrial<.15)
    qc=pd.DataFrame({'barcode':columns,'source_library':[PREFIXES[i] for i in assignment],
        'count_total':totals.astype(np.int64),'detected_genes':detected,
        'mitochondrial_fraction':mitochondrial,'passes_published_QC':paper_qc})
    qc.to_csv(out/'cell_QC.csv.gz',index=False)
    selectors={'all_exported':np.ones(n,dtype=bool),'published_QC':paper_qc}
    profiles={};captured=[]
    for mode,keep in selectors.items():
        number=np.bincount(assignment[keep],minlength=8)
        if (number==0).any():raise ValueError('QC excluded every cell in a library')
        for weighting in ['RNA_sum','equal_cell_mean']:
            profiles[(mode,weighting)]=[]
        for j,prefix in enumerate(PREFIXES):
            record=meta[GSMS[prefix]]
            if record.get('Sample_organism_ch1')!=['Mus musculus']:raise ValueError('Wrong source species')
            chars=record['Sample_characteristics_ch1']
            captured.append({'source_library':prefix,'gsm':GSMS[prefix],'QC':mode,
                'cells':int(number[j]),'original_title':record['Sample_title'][0],
                'original_characteristics':' | '.join(chars),
                'total_counts':int(totals[keep&(assignment==j)].sum())})
    with gzip.open(path,'rt') as stream:
        stream.readline()
        for row,line in enumerate(stream):
            gene,x=parse_count_row(line,n)
            if gene!=genes[row]:raise ValueError('Gene order changed between passes')
            for mode,keep in selectors.items():
                number=np.bincount(assignment[keep],minlength=8)
                profiles[(mode,'RNA_sum')].append(np.bincount(assignment[keep],weights=x[keep],minlength=8))
                profiles[(mode,'equal_cell_mean')].append(
                    np.bincount(assignment[keep],weights=x[keep]/totals[keep]*1e6,minlength=8)/number)
            if (row+1)%5000==0:print('PASS2',row+1,'genes',round(time.monotonic()-started,1),'seconds',flush=True)
    columns_meta=[];matrices=[]
    for (mode,weighting),rows in profiles.items():
        matrix=np.stack(rows);matrix=matrix/matrix.sum(axis=0)*1e6
        matrices.append(matrix)
        for prefix in PREFIXES:
            columns_meta.append({'profile':f'{prefix}:{mode}:{weighting}',
                'source_library':prefix,'QC':mode,'weighting':weighting,'gsm':GSMS[prefix]})
    frame=pd.DataFrame(np.column_stack(matrices),index=genes,columns=[r['profile'] for r in columns_meta])
    frame.index.name='gene';frame.to_csv(out/'library_profiles.csv.gz')
    pd.DataFrame(columns_meta).to_csv(out/'profile_evidence.csv',index=False)
    pd.DataFrame(captured).to_csv(out/'library_QC.csv',index=False)
    source_scripts='\n'.join((base/'author'/name).read_text() for name in ['R/Seurat_setup.R','R/Identify_Cell_Types_Manually.R'])
    audit={'cells_in_export':n,'genes_in_export':len(genes),'libraries':8,
        'cells_passing_published_QC':int(paper_qc.sum()),'profiles':frame.shape[1],
        'author_sample_table_rows':len(author),'author_sample_table_namespace':'ISO strict OOXML',
        'author_code_executed':False,'cell_type_truth_available':False,
        'QC_definition':'500 < detected genes < 8000; 2000 < counts < 125000; mitochondrial fraction < 0.15; source Seurat_setup.R bounds.',
        'author_Conditions_are_library_periods_not_cell_types':True,
        'label_gap':'Author scripts save labeled Seurat/RDS objects, but these objects/per-cell labels are absent from the inspected public repository tree and GEO attachment list.',
        'script_mentions_saved_labels':('StashIdent' in source_scripts and 'Cell.Types' in source_scripts),
        'limitations':['Captured/enriched libraries are not unbiased in-vivo proportions.',
            'Age and preparation/capture are confounded; no animal-level age p-values.',
            'Raw count totals and equal-cell means define different targets.',
            'QC sensitivity uses published thresholds; no independently known cell-type labels.'],
        'elapsed_seconds':round(time.monotonic()-started,2)}
    (out/'source_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':main()
