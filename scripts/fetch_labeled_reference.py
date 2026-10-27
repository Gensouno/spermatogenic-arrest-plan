#!/usr/bin/env python3
"""Bounded, cached E-MTAB-6946 source retrieval; never execute author code."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / 'data/raw/labeled_reference'
AUTHOR_COMMIT = '5d10fe6017651e4cc583d12ba0bac49c51717329'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''): h.update(block)
    return h.hexdigest()


def fetch(name, url, cap, download):
    path = CACHE/name
    if path.is_symlink(): raise ValueError('Symlink cache refused')
    if not path.exists():
        if not download: raise FileNotFoundError('Missing source; use --download: '+name)
        if not (os.environ.get('HTTPS_PROXY') and os.environ.get('https_proxy')):
            raise ValueError('Explicit HTTPS proxy required; no direct fallback')
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_name(path.name+'.download-part')
        if partial.is_symlink(): raise ValueError('Symlink partial refused')
        size = 0
        started=time.monotonic()
        request = urllib.request.Request(url, headers={'User-Agent': 'spermatogenic-arrest-plan-source-audit/1.0'})
        with urllib.request.urlopen(request, timeout=40) as response, partial.open('wb') as out:
            declared = response.headers.get('Content-Length')
            if declared and int(declared)>cap: raise ValueError('Declared file exceeds cap')
            for block in iter(lambda: response.read(1024**2), b''):
                size += len(block)
                if size>cap: raise ValueError('Stream exceeds cap')
                if time.monotonic()-started>360: raise TimeoutError('Six-minute transfer deadline; partial retained')
                out.write(block)
                if size%(100*1024**2)==0: print(name,size//1024**2,'MiB',round(time.monotonic()-started,1),'seconds',flush=True)
            if declared and size!=int(declared): raise ValueError('Incomplete download')
        partial.replace(path)
        print('Fetched', name, size, 'bytes', flush=True)
        time.sleep(.4)
    if path.stat().st_size>cap: raise ValueError('Cache exceeds cap')
    return {'path':str(path.relative_to(ROOT)), 'url':url,
            'bytes':path.stat().st_size, 'sha256':digest(path)}


def file_entries(value):
    if isinstance(value, dict):
        if 'files' in value:
            yield from file_entries_from_list(value['files'])
        for key, child in value.items():
            if key!='files': yield from file_entries(child)
    elif isinstance(value, list):
        for child in value: yield from file_entries(child)


def file_entries_from_list(value):
    if isinstance(value, list):
        for child in value: yield from file_entries_from_list(child)
    elif isinstance(value, dict) and 'path' in value:
        yield value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--download', action='store_true')
    args=p.parse_args()
    records=[fetch('study.json','https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-6946',2*1024**2,args.download),
             fetch('author_README.md',f'https://raw.githubusercontent.com/MarioniLab/Spermatogenesis2018/{AUTHOR_COMMIT}/README.md',1024**2,args.download)]
    study=json.loads((CACHE/'study.json').read_text())
    if study['accno']!='E-MTAB-6946': raise ValueError('Wrong source accession')
    entries=list(file_entries(study))
    print('Declared files:', json.dumps(entries,indent=2), flush=True)
    declared={r['path']:r for r in entries}
    for name in ['cell_metadata.txt','genes.tsv','E-MTAB-6946.sdrf.txt']:
        if name not in declared: raise ValueError('Required file not declared')
        row=fetch(name,'https://www.ebi.ac.uk/biostudies/files/E-MTAB-6946/'+name,15*1024**2,args.download)
        if row['bytes']!=declared[name]['size']: raise ValueError('Declared source size differs')
        records.append(row)
    out=ROOT/'results/labeled_reference';out.mkdir(parents=True,exist_ok=True)
    (out/'source_probe.json').write_text(json.dumps({'retrieved_utc':datetime.now(timezone.utc).isoformat(),
        'author_commit':AUTHOR_COMMIT,'assets':records,'declared_files':entries},indent=2)+'\n')


if __name__=='__main__': main()
