#!/usr/bin/env python3
"""Verify inputs and optional snapshot for the complete-P15 follow-up."""
import argparse
import json
from pathlib import Path
import re

from check_project import safe_path,verify_manifest
from emptydrops_range_client import digest

ROOT=Path(__file__).resolve().parents[1]


def input_rows(root):
    rows=json.loads((root/'data/emptydrops_p15_assets.json').read_text())['assets']
    rows+=[r for r in json.loads((root/'data/labeled_reference_assets.json').read_text())['assets']
        if Path(r['path']).name=='genes.tsv' or r.get('range',{}).get('library')=='do18195']
    if len(rows)!=8:raise ValueError('Source input inventory differs')
    return rows


def check_inputs(root):
    rows=input_rows(root)
    for row in rows:
        path=safe_path(root,row['path'])
        if path.stat().st_size!=row['bytes'] or digest(path)!=row['sha256']:raise ValueError('Source count/metadata cache differs: '+row['path'])
    fixed=[f'results/stage_reference/pooled_{m}_reference.csv.gz' for m in ['original_4','centroid_4','centroid_13']]
    # Reuse hashes fixed before this follow-up, rather than trusting new hashes.
    prior=json.loads((root/'labeled_reference_snapshot/design.json').read_text())['fixed_inputs_sha256']
    for relative in fixed:
        if digest(safe_path(root,relative))!=prior[relative]:raise ValueError('Original fixed reference differs')
    return len(rows)+len(fixed)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inputs',action='store_true');parser.add_argument('--inputs-only',action='store_true');args=parser.parse_args()
    if not args.inputs_only:
        print(verify_manifest(ROOT,'emptydrops_p15_snapshot'))
        for relative in ['docs/EMPTYDROPS_P15.md','emptydrops_p15_snapshot/README.md','README.md']:
            p=ROOT/relative
            for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
                link=link.strip('<>')
                if re.match(r'[A-Za-z][A-Za-z0-9+.-]*:',link) or link.startswith('#'):continue
                target=(p.parent/link.split('#')[0]).resolve()
                if not target.is_relative_to(ROOT) or not target.exists():raise ValueError('Broken link: '+relative)
    if args.inputs or args.inputs_only:print('Cached source/reference hashes checked:',check_inputs(ROOT))
    print('Integrity only; scientific validity and clean installation are not certified.')


if __name__=='__main__':main()
