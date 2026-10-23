#!/usr/bin/env python3
"""Download only explicit GEO-listed expansion assets under a size cap."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import urllib.request

from discover_expansion import parse_soft

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gse', required=True)
    parser.add_argument('--url', required=True)
    parser.add_argument('--max-mib', type=float, default=100)
    args = parser.parse_args()
    if not re.fullmatch(r'GSE[0-9]+', args.gse) or args.max_mib <= 0:
        raise ValueError('Valid GEO accession and positive size cap required')
    metadata = ROOT / f'data/raw/expansion_metadata/{args.gse}_family.soft.gz'
    series, samples = parse_soft(gzip.decompress(metadata.read_bytes()).decode())
    official = [u.replace('ftp://', 'https://') for record in [series, *samples.values()]
                for key, values in record.items() if 'supplementary_file' in key for u in values]
    if args.url not in official or not args.url.startswith('https://ftp.ncbi.nlm.nih.gov/geo/'):
        raise ValueError('Asset URL must exactly match this series official SOFT evidence')
    dest = ROOT / f'data/raw/expansion_assets/{args.gse}' / args.url.rsplit('/', 1)[-1]
    dest.parent.mkdir(parents=True, exist_ok=True)
    limit = int(args.max_mib * 1024**2)
    if not dest.exists():
        with urllib.request.urlopen(args.url, timeout=35) as response:
            length = response.headers.get('Content-Length')
            if length and int(length) > limit:
                raise ValueError(f'Asset exceeds configured size cap: {length} bytes')
            partial = dest.with_name(dest.name + '.download-part')
            if partial.is_symlink():
                raise ValueError('Refusing symlink partial')
            total = 0
            with partial.open('wb') as stream:
                for chunk in iter(lambda: response.read(1024 * 1024), b''):
                    total += len(chunk)
                    if total > limit:
                        raise ValueError('Asset exceeds configured size cap; partial retained')
                    stream.write(chunk)
                    if total % (10 * 1024**2) == 0:
                        print(dest.name, total // 1024**2, 'MiB', flush=True)
            if length and total != int(length):
                raise ValueError('Incomplete transfer; partial retained')
            partial.replace(dest)
    digest = hashlib.sha256()
    with dest.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    record = {'gse': args.gse, 'path': str(dest.relative_to(ROOT)), 'url': args.url,
              'bytes': dest.stat().st_size, 'sha256': digest.hexdigest(),
              'metadata_sha256': hashlib.sha256(metadata.read_bytes()).hexdigest()}
    manifest_path = ROOT / 'data/expansion_assets.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {'schema': 1, 'assets': []}
    existing = next((r for r in manifest['assets'] if r['path'] == record['path']), None)
    if existing and existing != record:
        raise ValueError('Existing expansion snapshot differs; do not silently replace its provenance')
    if not existing:
        manifest['assets'].append(record)
        manifest['updated_utc'] = datetime.now(timezone.utc).isoformat()
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    main()
