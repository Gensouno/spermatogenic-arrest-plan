#!/usr/bin/env python3
"""Collect fresh primary GEO evidence for optional expansion, not automatic inclusion."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def parse_soft(text):
    series, samples, target = {}, {}, None
    for line in text.splitlines():
        if line.startswith('^SERIES = '):
            target = series
        elif line.startswith('^SAMPLE = '):
            gsm = line.split(' = ', 1)[1]
            target = samples.setdefault(gsm, {})
        elif line.startswith('^'):
            target = None
        elif target is not None and line.startswith(('!Series_', '!Sample_')) and ' = ' in line:
            key, value = line.split(' = ', 1)
            target.setdefault(key.lstrip('!'), []).append(value)
    return series, samples


def acquire(gse):
    if not re.fullmatch(r'GSE[0-9]+', gse):
        raise ValueError('Invalid GEO series accession')
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/soft/{gse}_family.soft.gz'
    path = ROOT / f'data/raw/expansion_metadata/{gse}_family.soft.gz'
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        with urllib.request.urlopen(url, timeout=35) as response:
            data = response.read()
        text = gzip.decompress(data).decode('utf-8')
        if f'^SERIES = {gse}' not in text:
            raise ValueError('Unexpected SOFT content')
        path.write_bytes(data)
    data = path.read_bytes()
    series, samples = parse_soft(gzip.decompress(data).decode('utf-8'))
    for values in samples.values():
        values['Sample_supplementary_file'] = [url for key, urls in list(values.items())
            if key.startswith('Sample_supplementary_file') for url in urls]
    series_fields = ['title', 'summary', 'overall_design', 'type', 'pubmed_id', 'relation',
                     'supplementary_file', 'sample_id', 'status', 'last_update_date']
    sample_fields = ['title', 'source_name_ch1', 'organism_ch1', 'characteristics_ch1',
                     'library_strategy', 'library_source', 'library_selection',
                     'data_processing', 'supplementary_file', 'relation']
    return {'accession': gse, 'source_url': url, 'sha256': hashlib.sha256(data).hexdigest(),
            'series': {key: series.get('Series_' + key, []) for key in series_fields},
            'samples': {gsm: {key: values.get('Sample_' + key, []) for key in sample_fields}
                        for gsm, values in samples.items()},
            'inclusion': 'unreviewed; metadata discovery does not establish suitability'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--accessions', nargs='+', required=True)
    args = parser.parse_args()
    out = ROOT / 'results/data_expansion'
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for gse in dict.fromkeys(args.accessions):
        try:
            record = acquire(gse)
            records.append(record)
            (out / f'{gse}_metadata.json').write_text(json.dumps(record, indent=2) + '\n')
            print(gse, record['series']['title'], 'samples', len(record['samples']),
                  'series_files', record['series']['supplementary_file'], flush=True)
        except Exception as error:
            records.append({'accession': gse, 'error': str(error), 'inclusion': 'not evaluated'})
            print(gse, type(error).__name__, str(error), flush=True)
        time.sleep(.4)
    prior = {json.loads(path.read_text())['accession']: json.loads(path.read_text())
             for path in out.glob('GSE*_metadata.json')}
    prior.update({row['accession']: row for row in records})
    (out / 'discovery.json').write_text(json.dumps({'retrieved_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Targeted expansion search, not an exhaustive systematic review',
        'candidates': [prior[gse] for gse in sorted(prior)]}, indent=2) + '\n')


if __name__ == '__main__':
    main()
