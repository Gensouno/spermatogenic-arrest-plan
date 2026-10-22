#!/usr/bin/env python3
"""Frozen public-source assets: snapshot once, verify or download thereafter."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'data/analysis_assets.json'
COHORTS = ['GSE236388', 'GSE283904', 'GSE68283', 'GSE68284',
           'GSE101511', 'GSE294209', 'GSE38096', 'GSE42601']


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def definitions():
    rows = []
    def geo(gse, filename, local=None, kind='suppl'):
        prefix = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}'
        rows.append({'path': local or f'data/raw/{gse}/{filename}',
                     'url': f'{prefix}/{kind}/{filename}'})
    files = {'GSE236388': 'GSE236388_XAP5_counts.txt.gz',
             'GSE283904': 'GSE283904_Gene_count.txt.gz',
             'GSE68283': 'GSE68283_Raw_counts.txt.gz',
             'GSE68284': 'GSE68284_Raw_counts.txt.gz'}
    for gse in COHORTS + ['GSE112393']:
        filename = f'{gse}_family.soft.gz'
        geo(gse, filename, f'data/raw/metadata/{filename}', 'soft')
    for gse in COHORTS:
        geo(gse, files.get(gse, f'{gse}_RAW.tar'))
    for short in ['DGE', 'PerCellAttributes']:
        geo('GSE112393', f'GSE112393_MergedAdultMouseST25_{short}.txt.gz',
            f'data/raw/ref/GSE112393_{short}.txt.gz')
    filename = 'GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx'
    geo('GSE112393', filename, f'data/raw/ref/{filename}')
    filename = 'Mus_musculus.GRCm38.79.gtf.gz'
    rows.append({'path': f'data/ref/annotation/{filename}',
        'url': f'https://ftp.ensembl.org/pub/release-79/gtf/mus_musculus/{filename}'})
    filename = '12864_2016_2618_MOESM2_ESM.xlsx'
    rows.append({'path': f'data/raw/external_reference/{filename}',
        'url': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12864-016-2618-1/MediaObjects/' + filename})
    return rows


def safe_target(row):
    relative = Path(row['path'])
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError('Manifest destination must be a project-relative file')
    target = (ROOT / relative).resolve()
    if not target.is_relative_to(ROOT) or target == ROOT:
        raise ValueError('Asset path escapes project directory')
    if not row['url'].startswith('https://'):
        raise ValueError('Only explicit public HTTPS sources are supported')
    return target


def verify(row):
    target = safe_target(row)
    if not target.is_file():
        raise FileNotFoundError(row['path'])
    if target.stat().st_size != row['bytes'] or sha256(target) != row['sha256']:
        raise ValueError(f'Asset checksum mismatch: {row["path"]}; original retained, not overwritten')
    print('OK', row['path'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--snapshot', action='store_true', help='Create the frozen manifest from already acquired sources; does not overwrite it')
    action.add_argument('--verify', action='store_true')
    action.add_argument('--download', action='store_true', help='Download only missing assets using resumable curl, then verify SHA-256')
    action.add_argument('--import-from', type=Path, help='Copy frozen assets from another project directory; verify before copying, never replace existing input')
    args = parser.parse_args()
    if args.snapshot:
        if MANIFEST.exists():
            raise FileExistsError('Frozen manifest already exists; review explicitly rather than silently replacing hashes')
        rows = definitions()
        for row in rows:
            path = safe_target(row)
            row.update(bytes=path.stat().st_size, sha256=sha256(path))
            if '/geo/series/' in row['url'] and '/suppl/' in row['url']:
                gse = row['url'].split('/')[6]
                record = gzip.decompress((ROOT / f'data/raw/metadata/{gse}_family.soft.gz').read_bytes()).decode()
                if row['url'].replace('https://', 'ftp://') not in record:
                    raise ValueError(f'Supplement URL not listed in official series SOFT: {row["url"]}')
        MANIFEST.write_text(json.dumps({'schema': 1, 'scope': 'Eight selected cohorts, atlas and two reference supplements; no unused discovery data', 'assets': rows}, indent=2) + '\n')
        print('Wrote', MANIFEST, 'assets:', len(rows))
        return
    manifest = json.loads(MANIFEST.read_text())
    for row in manifest['assets']:
        target = safe_target(row)
        if args.import_from is not None and not target.exists():
            source = (args.import_from / row['path']).resolve()
            if source.stat().st_size != row['bytes'] or sha256(source) != row['sha256']:
                raise ValueError(f'Local source differs from frozen input: {source}')
            target.parent.mkdir(parents=True, exist_ok=True)
            # Exclusive creation protects existing or symlink partial files.
            partial = target.with_name(target.name + '.import-part')
            with source.open('rb') as src, partial.open('xb') as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
            if partial.stat().st_size != row['bytes'] or sha256(partial) != row['sha256']:
                raise ValueError(f'Local copy checksum mismatch: {partial}')
            partial.replace(target)
        if args.download and not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            partial = target.with_name(target.name + '.download-part')
            if partial.is_symlink():
                raise ValueError('Refusing symlink download partial')
            # curl respects the caller's proxy environment; no local proxy is
            # hardcoded into the published reproducibility instructions.
            subprocess.run(['curl', '--fail', '--location', '--retry', '3',
                '--continue-at', '-', '--connect-timeout', '30',
                '--output', str(partial), row['url']], check=True)
            if partial.stat().st_size != row['bytes'] or sha256(partial) != row['sha256']:
                raise ValueError(f'Download differs from frozen source: {partial}; not promoted to input')
            partial.replace(target)
        verify(row)


if __name__ == '__main__':
    main()
