#!/usr/bin/env python3
"""New-source provenance and snapshot integrity only."""
import argparse
import json
from pathlib import Path
from check_project import safe_path,verify_manifest
from freeze_replicability_inputs import digest
ROOT=Path(__file__).resolve().parents[1]


def check_inputs(root):
    rows=json.loads((root/'data/new_model_inputs.json').read_text())['assets']
    if len({r['path'] for r in rows})!=len(rows):raise ValueError('Duplicate input paths')
    for row in rows:
        p=safe_path(root,row['path'])
        if p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Frozen new-source input differs: '+row['path'])
    return len(rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',action='store_true');p.add_argument('--inputs-only',action='store_true');args=p.parse_args()
    if not args.inputs_only:print(verify_manifest(ROOT,'new_models_snapshot'))
    if args.inputs or args.inputs_only:print('Frozen new-source inputs:',check_inputs(ROOT))
    print('Integrity only; source animal independence and scientific validity not established.')


if __name__=='__main__':main()
