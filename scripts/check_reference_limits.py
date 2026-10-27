#!/usr/bin/env python3
"""Check complete reference-limit source/prepared input and compact result integrity."""
import argparse
import json
from pathlib import Path
import re

from check_project import safe_path,verify_manifest
from freeze_reference_limits_inputs import digest

ROOT=Path(__file__).resolve().parents[1]


def check_inputs(root):
    rows=json.loads((root/'data/reference_limits_inputs.json').read_text())['assets']
    rows+=json.loads((root/'data/reference_limits_sorted_inputs.json').read_text())['assets']
    for row in rows:
        path=safe_path(root,row['path'])
        if path.stat().st_size!=row['bytes'] or digest(path)!=row['sha256']:raise ValueError('Frozen source/prepared input differs: '+row['path'])
    hashes={r['path']:r['sha256'] for r in rows}
    for name,keys in [('technical_design.json',['fixed_inputs_sha256']),
        ('transfer_design.json',['fixed_prepared_inputs_sha256','generated_candidate_reference_sha256']),
        ('legacy_sorted_design.json',['fixed_candidate_reference_sha256'])]:
        path=root/'reference_limits_snapshot'/name
        if not path.exists():continue
        document=json.loads(path.read_text())
        for key in keys:
            for relative,checksum in document.get(key,{}).items():
                if relative in hashes and hashes[relative]!=checksum:raise ValueError('Conflicting derived/source fingerprints')
                if digest(safe_path(root,relative))!=checksum:raise ValueError('Prepared/candidate input differs: '+relative)
                hashes[relative]=checksum
    return len(hashes)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inputs',action='store_true');parser.add_argument('--inputs-only',action='store_true');args=parser.parse_args()
    if not args.inputs_only:
        print(verify_manifest(ROOT,'reference_limits_snapshot'))
        for relative in ['README.md','docs/REFERENCE_LIMITS.md','reference_limits_snapshot/README.md']:
            p=ROOT/relative
            for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
                link=link.strip('<>')
                if re.match(r'[A-Za-z][A-Za-z0-9+.-]*:',link) or link.startswith('#'):continue
                target=(p.parent/link.split('#')[0]).resolve()
                if not target.is_relative_to(ROOT) or not target.exists():raise ValueError('Broken link '+relative)
    if args.inputs or args.inputs_only:print('Frozen raw/prepared source hashes checked:',check_inputs(ROOT))
    print('Integrity only; not scientific validity or clean-installation certification.')


if __name__=='__main__':main()
