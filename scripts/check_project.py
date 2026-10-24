#!/usr/bin/env python3
"""Offline integrity checks for all available project snapshots, including local work."""
import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
SNAPSHOTS=['results_snapshot','expansion_snapshot','followup_snapshot','reference_audit_snapshot']


def safe_path(root,relative):
    p=Path(relative)
    if p.is_absolute() or '..' in p.parts or not p.parts:raise ValueError('Unsafe manifest path')
    target=root/p
    if target.is_symlink() or not target.resolve().is_relative_to(root.resolve()):raise ValueError('Manifest path escape/symlink')
    return target


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(root,snapshot):
    base=safe_path(root,snapshot);manifest=json.loads((base/'manifest.json').read_text())
    names=[r['path'] for r in manifest['files']]
    actual={str(p.relative_to(base)) for p in base.rglob('*') if p.is_file()
            and str(p.relative_to(base)) not in ['manifest.json','.DS_Store']}
    if len(names)!=len(set(names)) or set(names)!=actual:raise ValueError('Snapshot inventory differs')
    for row in manifest['files']:
        p=safe_path(base,row['path'])
        if p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Snapshot content changed: '+row['path'])
    specs=dict(manifest.get('source_sha256',{}));specs.update(manifest.get('input_specification_sha256',{}))
    conventional={'results_snapshot':[('data/analysis_assets.json','asset_manifest_sha256'),('data/samples.tsv','sample_sheet_sha256')],
        'expansion_snapshot':[('data/expansion_assets.json','asset_manifest_sha256'),('data/expansion_review.json','review_sha256')]}
    specs.update({p:manifest[k] for p,k in conventional.get(snapshot,[])})
    for relative,checksum in specs.items():
        if digest(safe_path(root,relative))!=checksum:raise ValueError('Source/specification changed: '+relative)
    return {'snapshot':snapshot,'files':len(names),'source_and_specification_checks':len(specs),'status':'verified'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',action='store_true',help='Also hash locally cached raw/prepared/reference inputs; may be slow')
    args=parser.parse_args()
    reports=[verify_manifest(ROOT,n) for n in SNAPSHOTS if (ROOT/n/'manifest.json').exists()]
    for r in reports:print(r)
    documents=['README.md','docs/DATA_EXPANSION.md','docs/FOLLOWUP_LRH1.md','docs/REFERENCE_TRANSFER.md']
    documents += [n+'/README.md' for n in SNAPSHOTS]
    for relative in documents:
        p=ROOT/relative
        if not p.exists():continue
        for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            link=link.strip('<>')
            if re.match(r'[A-Za-z][A-Za-z0-9+.-]*:',link) or link.startswith('#'):continue
            target=(p.parent/link.split('#')[0]).resolve()
            if not target.is_relative_to(ROOT) or not target.exists():raise ValueError('Broken document link: '+relative)
    raw_checks=0
    if args.inputs:
        for name in ['data/followup_assets.json','data/age_reference_assets.json']:
            if not (ROOT/name).exists():continue
            for row in json.loads((ROOT/name).read_text())['assets']:
                p=safe_path(ROOT,row['path'])
                if p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Raw cache differs: '+row['path'])
                raw_checks+=1
        p=ROOT/'reference_audit_snapshot/design.json'
        if p.exists():
            for name,checksum in json.loads(p.read_text())['fixed_inputs_sha256'].items():
                if digest(safe_path(ROOT,name))!=checksum:raise ValueError('Prepared/reference input changed: '+name)
                raw_checks+=1
    print('Offline project checks passed; cached input hashes checked:',raw_checks)
    print('This is file integrity, not a scientific validity or clean-installation certificate.')


if __name__=='__main__':main()
