#!/usr/bin/env python3
"""Read-only Git preparation and publication guards; never commit or push."""
import argparse
import fnmatch
import json
from pathlib import Path
import re
import subprocess
from paper_support import ROOT,safe_path

PRIVATE=['paper/*','research/*','.env','.env.*','.download-proxy.local','.github-sync.local.json',
    '.paper-tools/*','.paper-private/*','data/raw/*','data/bulk/*','data/interim/*','logs/*','.venv*','.Rlib*']
SECRET=[re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
        re.compile(rb'gh[pousr]_[A-Za-z0-9]{30,}'),re.compile(rb'github_pat_[A-Za-z0-9_]{40,}'),
        re.compile(rb'AKIA[0-9A-Z]{16}'),re.compile(rb'sk-[A-Za-z0-9_-]{32,}')]


def git(root,*args):
    return subprocess.run(['git','-C',str(root),*args],capture_output=True,check=True).stdout


def private_path(path):return any(fnmatch.fnmatchcase(path,p) for p in PRIVATE)


def inspect(root):
    root=Path(root).resolve();issues=[]
    files=git(root,'ls-files','--cached','--others','--exclude-standard','-z').decode().split('\0')
    tracked=set(git(root,'ls-files','-z').decode().split('\0'))
    changed=git(root,'status','--porcelain=v1','-z').decode().split('\0')
    for name in sorted(set(files)-{''}):
        if private_path(name):issues.append({'level':'ERROR','code':'PROTECTED_PATH','path':name});continue
        try:
            p=safe_path(root,name)
            if p.is_symlink() or not p.is_file():continue
            if p.stat().st_size>20*1024*1024:issues.append({'level':'WARNING','code':'LARGE_FILE','path':name})
            if p.stat().st_size<=5*1024*1024:
                content=p.read_bytes()
                if b'\0' not in content[:8192] and any(pattern.search(content) for pattern in SECRET):
                    issues.append({'level':'ERROR','code':'POSSIBLE_SECRET','path':name})
        except (OSError,ValueError):issues.append({'level':'WARNING','code':'PATH_REVIEW','path':name})
    return {'operation':'read_only_preparation','HEAD':git(root,'rev-parse','HEAD').decode().strip(),
        'branch':git(root,'branch','--show-current').decode().strip(),
        'remote':git(root,'remote','get-url','origin').decode().strip() if git(root,'remote').strip() else None,
        'changed_status_records':[x for x in changed if x],
        'private_material_tracked':[p for p in tracked if p and private_path(p)],
        'issues':issues,'commit_suggestion':'feat: add evidence-linked manuscript synchronization checks',
        'authorization':'No commit, push, publication or visibility change performed. Verify actual user authorization and release scope before any mutation.',
        'limits':'Pattern scanning is not a complete secret, privacy, licensing or unpublished-content audit. Ignored local manuscript changes are intentionally not proposed for the public remote.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--report',type=Path);args=p.parse_args()
    try:report=inspect(args.root)
    except (OSError,subprocess.CalledProcessError,ValueError) as exc:print('ERROR GIT_CHECK:',type(exc).__name__);return 1
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(report,indent=2)+'\n')
    for i in report['issues']:print(i['level'],i['code'],i['path'])
    print('INFO',len(report['changed_status_records']),'Git status records; no commit/push performed.')
    print('INFO Unpublished paper/research stay local; manual release review and scoped authorization required.')
    return int(any(i['level']=='ERROR' for i in report['issues']))


if __name__=='__main__':raise SystemExit(main())
