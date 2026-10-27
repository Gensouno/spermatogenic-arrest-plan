#!/usr/bin/env python3
"""Verify the replicability snapshot and optional frozen local inputs."""
import argparse
import json
from pathlib import Path
from check_project import verify_manifest,safe_path
from freeze_replicability_inputs import digest
ROOT=Path(__file__).resolve().parents[1]


def check_inputs(root):
    assets=json.loads((root/'data/replicability_inputs.json').read_text())['assets']
    if len({r['path'] for r in assets})!=len(assets):raise ValueError('Duplicate frozen path')
    for row in assets:
        p=safe_path(root,row['path'])
        if p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Frozen input differs: '+row['path'])
    return len(assets)


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',action='store_true');p.add_argument('--inputs-only',action='store_true');args=p.parse_args()
    if not args.inputs_only:print(verify_manifest(ROOT,'replicability_snapshot'))
    if args.inputs or args.inputs_only:print('Frozen inputs checked:',check_inputs(ROOT))
    print('File integrity only; no scientific validity claim.')


if __name__=='__main__':main()
