#!/usr/bin/env python3
"""Bounded GEO follow-up retrieval, isolated from both published analyses."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import time
import urllib.request

from discover_expansion import parse_soft

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'data/followup_assets.json'
FILES = {
    'GSE180204': ['GSE180204_Agrimson_Lrh1_Supplemental_Tables.xlsx',
                  'GSE180204_Lrh1_p1_p3_RnaSeq_Agrimson.rdata.gz'],
    'GSE181132': ['GSE181132_RAW.tar'],
}
LIMIT = 20 * 1024**2


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def get_file(url, path, expected=None):
    if path.is_symlink():
        raise ValueError('Refusing a symlink input')
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_name(path.name + '.download-part')
        if partial.is_symlink():
            raise ValueError('Refusing a symlink partial')
        total = 0
        with urllib.request.urlopen(url, timeout=35) as response, partial.open('wb') as stream:
            length = response.headers.get('Content-Length')
            if length and int(length) > LIMIT:
                raise ValueError('Input exceeds 20 MiB follow-up cap')
            for chunk in iter(lambda: response.read(1024**2), b''):
                total += len(chunk)
                if total > LIMIT:
                    raise ValueError('Input exceeds 20 MiB follow-up cap')
                stream.write(chunk)
        if length and total != int(length):
            raise ValueError('Incomplete transfer; partial retained')
        if expected and (total != expected['bytes'] or digest(partial) != expected['sha256']):
            raise ValueError('Downloaded input differs from frozen source; partial retained')
        partial.replace(path)
        time.sleep(.4)
    record = {'path': str(path.relative_to(ROOT)), 'url': url,
              'bytes': path.stat().st_size, 'sha256': digest(path)}
    if expected and any(record[k] != expected[k] for k in record):
        raise ValueError('Cached input differs from frozen source')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true', help='Explicitly fetch missing inputs; otherwise require cache')
    args = parser.parse_args()
    prior = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else None
    expected = {r['path']: r for r in prior['assets']} if prior else {}
    records, metadata = [], {}
    for gse, names in FILES.items():
        base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}'
        relative = f'data/raw/followup/metadata/{gse}_family.soft.gz'
        requests = [(base + f'/soft/{gse}_family.soft.gz', relative)]
        requests += [(base + '/suppl/' + name, f'data/raw/followup/{gse}/{name}') for name in names]
        for index, (url, relative) in enumerate(requests):
            path = ROOT / relative
            if not args.download and not path.exists():
                raise FileNotFoundError('Missing cache; use --download explicitly: ' + relative)
            if index:
                official = metadata[gse]['series'].get('Series_supplementary_file', [])
                if url.replace('https://', 'ftp://') not in official:
                    raise ValueError('Asset is not listed in current official SOFT')
            row = get_file(url, path, expected.get(relative))
            row['gse'] = gse
            records.append(row)
            if index == 0:
                text = gzip.decompress(path.read_bytes()).decode()
                if f'^SERIES = {gse}' not in text:
                    raise ValueError('Unexpected series metadata')
                series, samples = parse_soft(text)
                metadata[gse] = {'series': series, 'samples': samples,
                                 'url': url, 'sha256': row['sha256']}
            print(gse, path.name, row['bytes'], 'bytes', flush=True)
    if prior and prior['assets'] != records:
        raise ValueError('Frozen follow-up inventory changed; review rather than overwrite')
    if not prior:
        MANIFEST.write_text(json.dumps({'schema': 1,
            'retrieved_utc': datetime.now(timezone.utc).isoformat(),
            'scope': 'Targeted follow-up of GSE180204 and GSE181132, not exhaustive search',
            'assets': records}, indent=2) + '\n')
    out = ROOT / 'results/followup'; out.mkdir(parents=True, exist_ok=True)
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('Reconciled:', len(FILES), 'series;', len(records), 'files;', sum(r['bytes'] for r in records), 'bytes')


if __name__ == '__main__':
    main()
