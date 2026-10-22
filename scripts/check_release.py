#!/usr/bin/env python3
"""Pre-push checks of the reviewed Git index, documentation and result snapshot."""
from pathlib import Path
import re
import subprocess

from export_snapshot import verify_snapshot

ROOT = Path(__file__).resolve().parents[1]
BLOCKED_PREFIXES = ['.venv', '.Rlib', '.cache', 'data/raw', 'data/ref', 'data/bulk',
                    'data/interim', 'results/', 'logs/', '__pycache__']
SECRET = re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{20,}|sk-[A-Za-z0-9]{30,}')
DOCS = ['README.md', 'data/README.md', 'results_snapshot/README.md',
        'docs/METHODS.md', 'docs/RESULTS.md', 'docs/REPRODUCIBILITY.md',
        'docs/HISTORY.md', 'docs/INTERVIEW.md', 'notes/README.md']


def main():
    subprocess.run(['git', 'diff', '--quiet'], cwd=ROOT, check=True)
    verify_snapshot()
    staged = subprocess.check_output(['git', 'ls-files', '-s', '-z'], cwd=ROOT).decode().split('\0')
    paths = set()
    for entry in filter(None, staged):
        metadata, relative = entry.split('\t', 1)
        mode = metadata.split()[0]
        paths.add(relative)
        if mode == '120000':
            raise ValueError(f'Do not publish local environment symlinks: {relative}')
        if any(relative == prefix.rstrip('/') or relative.startswith(prefix.rstrip('/') + '/') for prefix in BLOCKED_PREFIXES):
            raise ValueError(f'Generated input/environment accidentally tracked: {relative}')
        path = ROOT / relative
        if path.stat().st_size > 2 * 1024 * 1024:
            raise ValueError(f'Review large tracked file explicitly: {relative}')
        if path.suffix in ['.py', '.R', '.sh', '.md', '.json', '.csv', '.tsv', '.txt']:
            text = path.read_text()
            if SECRET.search(text):
                raise ValueError(f'Potential credential material: {relative}; contents withheld')
            if relative.startswith('results_snapshot/') and re.search(r'/Users/|/private/(?:var|tmp)/|/var/folders/', text):
                raise ValueError(f'Local machine path in public snapshot: {relative}')
    for relative in DOCS:
        for match in re.finditer(r'!?\[[^\]]*\]\(([^)]+)\)', (ROOT / relative).read_text()):
            target = match.group(1).strip('<>')
            if re.match(r'[a-zA-Z][a-zA-Z0-9+.-]*:', target) or target.startswith('#'):
                continue
            target = target.split('#', 1)[0]
            resolved = (ROOT / relative).parent.joinpath(target).resolve()
            if not resolved.is_relative_to(ROOT) or not resolved.exists():
                raise ValueError(f'Broken active documentation link: {relative}: {target}')
            if str(resolved.relative_to(ROOT)) not in paths:
                raise ValueError(f'Active documentation target is not in Git index: {target}')
    print('Release index checks OK:', len(paths), 'tracked files; no raw data, dependency links, large matrices or detected credential patterns')


if __name__ == '__main__':
    main()
