#!/usr/bin/env python3
"""Fetch official GEO SOFT records and join selected samples to source evidence."""
import gzip
import hashlib
import json
import re
from pathlib import Path
import urllib.request

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    dest = ROOT / 'data/raw/metadata'
    dest.mkdir(parents=True, exist_ok=True)
    all_rows, sources = [], []
    for gse in ['GSE236388','GSE283904','GSE68283','GSE68284','GSE101511','GSE294209','GSE38096','GSE42601']:
        url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/soft/{gse}_family.soft.gz'
        path = dest / f'{gse}_family.soft.gz'
        if not path.exists():
            with urllib.request.urlopen(url, timeout=20) as response:
                content = response.read()
            text = gzip.decompress(content).decode()
            if f'^SERIES = {gse}' not in text:
                raise ValueError(f'Incorrect accession content {gse}')
            path.write_bytes(content)
        text = gzip.decompress(path.read_bytes()).decode()
        records, current = {}, None
        for line in text.splitlines():
            if line.startswith('^SAMPLE = '):
                current = line.split(' = ',1)[1]
                records[current] = {}
            elif line.startswith('^'):
                current = None
            elif current and line.startswith(('!Sample_title = ', '!Sample_source_name_ch1 = ', '!Sample_characteristics_ch1 = ')):
                key, value = line.split(' = ',1)
                records[current].setdefault(key,[]).append(value)
        sheet_path = ROOT / (f'data/interim/arrays/{gse}/samples.tsv' if gse in ['GSE38096','GSE42601']
                             else f'data/bulk/{gse}_samples.tsv')
        selected = pd.read_csv(sheet_path,sep='\t')
        for _, sample in selected.iterrows():
            record = records[sample.sample_id]
            title = ' | '.join(record['!Sample_title'])
            normalize = lambda value: re.sub(r'\s*,\s*', ',', ' '.join(value.split()))
            curated_title = sample.title
            if gse == 'GSE101511':
                # Earlier inventory prepended genotype to the actual GEO title.
                # Verify the title suffix plus genotype independently below.
                curated_title = re.sub(r'^(Hnrnpgt-/-|wild type)\s+', '', curated_title)
                genotype = ' | '.join(record.get('!Sample_characteristics_ch1', []))
                expected = 'genotype: Hnrnpgt-/-' if sample.group == 'KO' else 'genotype: wild type'
                if expected not in genotype:
                    raise ValueError(f'Genotype mismatch {sample.sample_id}')
            if normalize(title) != normalize(curated_title):
                raise ValueError(f'Title mismatch {sample.sample_id}: {title} vs {sample.title}')
            all_rows.append({'gse':gse,'gsm':sample.sample_id,'matrix_sample':sample.matrix_sample,
                'assigned_group':sample.group,'official_title':title,
                'official_source':' | '.join(record.get('!Sample_source_name_ch1',[])),
                'official_characteristics':' | '.join(record.get('!Sample_characteristics_ch1',[])),
                'title_matches_after_whitespace_normalization':True,'url':url})
        sources.append({'gse':gse,'url':url,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                        'selected_samples':len(selected),'series_samples':len(records)})
        print(gse, 'verified selected sample IDs and titles:', len(selected), flush=True)
    out = ROOT / 'results/import_qc'
    pd.DataFrame(all_rows).to_csv(out / 'official_sample_evidence.csv',index=False)
    (out / 'metadata_sources.json').write_text(json.dumps(sources,indent=2)+'\n')


if __name__ == '__main__':
    main()
