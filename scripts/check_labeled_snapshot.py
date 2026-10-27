#!/usr/bin/env python3
"""Verify optional labelled-library snapshot without changing old snapshot contracts."""
import argparse
import json
from pathlib import Path
import re

from check_project import digest,safe_path,verify_manifest

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',action='store_true');a=p.parse_args()
    print(verify_manifest(ROOT,'labeled_reference_snapshot'))
    for relative in ['README.md','docs/LABELED_REFERENCE.md','labeled_reference_snapshot/README.md']:
        document=ROOT/relative
        for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',document.read_text()):
            link=link.strip('<>')
            if re.match(r'[A-Za-z][A-Za-z0-9+.-]*:',link) or link.startswith('#'):continue
            target=(document.parent/link.split('#')[0]).resolve()
            if not target.is_relative_to(ROOT) or not target.exists():raise ValueError('Broken documentation link: '+relative)
    if a.inputs:
        rows=json.loads((ROOT/'data/labeled_reference_assets.json').read_text())['assets']
        rows+=json.loads((ROOT/'data/emptydrops_label_assets.json').read_text())['assets']
        for row in rows:
            path=safe_path(ROOT,row['path'])
            if path.stat().st_size!=row['bytes'] or digest(path)!=row['sha256']:raise ValueError('Source cache differs')
        fixed=json.loads((ROOT/'labeled_reference_snapshot/design.json').read_text())['fixed_inputs_sha256']
        for relative,checksum in fixed.items():
            if digest(safe_path(ROOT,relative))!=checksum:raise ValueError('Frozen reference differs')
        print('Cached source/reference input hashes checked:',len(rows)+len(fixed))
    print('File integrity only; not a scientific validity or clean-installation certificate.')


if __name__=='__main__':main()
