#!/usr/bin/env python3
"""Verify independent-source snapshot without changing prior snapshot contracts."""
import argparse
import json
from pathlib import Path
import re

from check_project import safe_path,verify_manifest
from export_leptotene import digest

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inputs',action='store_true');args=parser.parse_args()
    print(verify_manifest(ROOT,'leptotene_snapshot'))
    for relative in ['README.md','docs/LEPTOTENE_FOLLOWUP.md','leptotene_snapshot/README.md']:
        p=ROOT/relative
        for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            link=link.strip('<>')
            if re.match(r'[A-Za-z][A-Za-z0-9+.-]*:',link) or link.startswith('#'):continue
            target=(p.parent/link.split('#')[0]).resolve()
            if not target.is_relative_to(ROOT) or not target.exists():raise ValueError('Broken documentation link: '+relative)
    if args.inputs:
        rows=json.loads((ROOT/'data/leptotene_assets.json').read_text())['assets']
        for row in rows:
            p=safe_path(ROOT,row['path'])
            if p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Source asset differs')
        fixed=json.loads((ROOT/'leptotene_snapshot/design.json').read_text())['fixed_inputs_sha256']
        for name,checksum in fixed.items():
            if digest(safe_path(ROOT,name))!=checksum:raise ValueError('Frozen reference/mapping/paper input differs')
        print('Cached source/reference/mapping/paper hashes checked:',len(rows)+len(fixed))
    print('Integrity only; not a scientific validity or clean-installation certificate.')


if __name__=='__main__':main()
