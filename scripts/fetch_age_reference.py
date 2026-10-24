#!/usr/bin/env python3
"""Bounded primary-source age-reference retrieval; no author code execution."""
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
COMMIT = '649d8aa4976915c834f195b79eff161cd10c38c6'
BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121904'
MATRIX = 'GSE121904_testis_maturation_scRNAseq_raw_mapped_counts.csv.gz'
AUTHOR = ['README.md', 'LICENSE.md', 'doc/20180530_sample_info.xlsx',
          'R/Seurat_setup.R', 'R/Identify_Cell_Types_Manually.R']


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''): h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    manifest_path = ROOT/'data/age_reference_assets.json'
    frozen = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
    expected = {r['path']: r for r in frozen['assets']} if frozen else {}
    requests = [('GSE121904_family.soft.gz', BASE+'/soft/GSE121904_family.soft.gz', 1024**2),
                (MATRIX, BASE+'/suppl/'+MATRIX, 80*1024**2)]
    requests += [('author/'+name, f'https://raw.githubusercontent.com/nyuhuyang/scRNAseq-SSCs/{COMMIT}/{name}',
                  1024**2) for name in AUTHOR]
    records = []
    series = None
    for relative, url, cap in requests:
        path = ROOT/'data/raw/age_reference'/relative
        key = str(path.relative_to(ROOT))
        if path.is_symlink(): raise ValueError('Refusing symlink cache')
        if relative == MATRIX and url.replace('https://','ftp://') not in series['Series_supplementary_file']:
            raise ValueError('Matrix URL not present in official metadata')
        if not path.exists():
            if not args.download: raise FileNotFoundError('Missing cache; use --download: '+key)
            path.parent.mkdir(parents=True, exist_ok=True)
            partial = path.with_name(path.name+'.download-part')
            if partial.is_symlink(): raise ValueError('Refusing symlink partial')
            start=time.monotonic();total=0
            with urllib.request.urlopen(url,timeout=35) as response, partial.open('wb') as stream:
                length=response.headers.get('Content-Length')
                if length and int(length)>cap: raise ValueError('Exceeds size cap')
                for chunk in iter(lambda: response.read(1024**2),b''):
                    total+=len(chunk)
                    if total>cap: raise ValueError('Exceeds size cap')
                    stream.write(chunk)
                    if total%(10*1024**2)==0: print(relative,total//1024**2,'MiB',flush=True)
                if length and total!=int(length): raise ValueError('Incomplete transfer')
            if key in expected and (total!=expected[key]['bytes'] or sha256(partial)!=expected[key]['sha256']):
                raise ValueError('Input differs from frozen source; partial retained')
            partial.replace(path);time.sleep(.4)
            print('Fetched',relative,round(time.monotonic()-start,1),'seconds',flush=True)
        row={'path':key,'url':url,'bytes':path.stat().st_size,'sha256':sha256(path)}
        if key in expected and row!=expected[key]: raise ValueError('Frozen source changed')
        records.append(row)
        if relative.endswith('_family.soft.gz'):
            text=gzip.decompress(path.read_bytes()).decode()
            if '^SERIES = GSE121904' not in text: raise ValueError('Wrong GEO series')
            series,samples=parse_soft(text)
            if len(samples)!=8: raise ValueError('Expected eight library records')
    if frozen and records!=frozen['assets']: raise ValueError('Frozen inventory changed')
    if not frozen:
        manifest_path.write_text(json.dumps({'schema':1,'retrieved_utc':datetime.now(timezone.utc).isoformat(),
            'scope':'One mouse developmental count matrix, eight GEO sample records and five author-source files; no per-cell truth labels assumed.',
            'author_commit':COMMIT,'author_code_executed':False,'assets':records},indent=2)+'\n')
    out=ROOT/'results/age_reference';out.mkdir(parents=True,exist_ok=True)
    (out/'metadata.json').write_text(json.dumps({'series':series,'samples':samples},indent=2)+'\n')
    print('Verified',len(records),'files;',sum(r['bytes'] for r in records),'bytes')


if __name__=='__main__': main()
