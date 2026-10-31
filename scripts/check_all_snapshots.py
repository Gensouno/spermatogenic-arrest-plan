#!/usr/bin/env python3
"""Discover every root snapshot and check its immutable artifacts and bound sources."""
import argparse
import json
from pathlib import Path
from check_project import verify_manifest


def check(root):
    snapshots=sorted(p for p in root.glob('*_snapshot') if p.is_dir())
    if not snapshots:raise ValueError('No snapshot directories found')
    rows=[verify_manifest(root,p.name) for p in snapshots]
    return dict(scope='Artifact/source integrity only, not numerical reruns or biological validation.',
        snapshots=len(rows),files=sum(r['files'] for r in rows),
        source_specification_checks=sum(r['source_and_specification_checks'] for r in rows),details=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--report',type=Path);a=parser.parse_args()
    result=check(Path(__file__).resolve().parents[1])
    if a.report:
        a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='details'},indent=2))
