"""Small strict YAML, path, digest and structured-result helpers for paper tools."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
LOCAL=ROOT/'.paper-tools/packages'
if LOCAL.is_dir():sys.path.insert(0,str(LOCAL))
try:
    import yaml
except ImportError as exc:
    raise SystemExit('Missing PyYAML. Install envs/requirements-paper-check.txt or use envs/paper_python.sh.') from exc


class StrictLoader(yaml.SafeLoader):pass


def unique_mapping(loader,node,deep=False):
    result={}
    for k,v in node.value:
        key=loader.construct_object(k,deep=deep)
        if key in result:raise ValueError('Duplicate YAML key: '+str(key))
        result[key]=loader.construct_object(v,deep=deep)
    return result


StrictLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,unique_mapping)


def read_yaml(path):
    try:value=yaml.load(path.read_text(encoding='utf-8'),Loader=StrictLoader)
    except yaml.YAMLError as exc:raise ValueError('Invalid YAML syntax in '+path.name) from exc
    if not isinstance(value,dict):raise ValueError('YAML document must be a mapping')
    return value


def safe_path(root,relative):
    p=Path(relative)
    if p.is_absolute() or '..' in p.parts or not p.parts:raise ValueError('Unsafe relative path')
    target=root/p
    if not target.resolve().is_relative_to(root.resolve()) or target.is_symlink():raise ValueError('Path escape/symlink')
    return target


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def frontmatter(path):
    text=path.read_text(encoding='utf-8')
    if not text.startswith('---\n'):raise ValueError('Chapter requires YAML frontmatter')
    parts=text.split('\n---\n',1)
    if len(parts)!=2:raise ValueError('Unclosed chapter frontmatter')
    try:meta=yaml.load(parts[0][4:],Loader=StrictLoader)
    except yaml.YAMLError as exc:raise ValueError('Invalid frontmatter YAML in '+path.name) from exc
    if not isinstance(meta,dict):raise ValueError('Chapter frontmatter must be a mapping')
    return meta,parts[1]


def result_value(root,spec):
    path=safe_path(root,spec['source'])
    if spec['format']=='json':
        obj=json.loads(path.read_text())
        for key in spec['pointer']:obj=obj[key]
        return obj
    if spec['format']!='csv':raise ValueError('Unsupported numeric source format')
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
    where=spec.get('where',{})
    rows=[r for r in rows if all(str(r.get(k))==str(v) for k,v in where.items())]
    op=spec.get('operation','value')
    if op=='row_count':return len(rows)
    if op=='unique_count':return len({r[spec['column']] for r in rows})
    if op=='sum':return sum(float(r[spec['column']]) for r in rows)
    if op=='sum_columns':return sum(float(r[k]) for r in rows for k in spec['columns'])
    if op in ['min','max']:
        if not rows:raise ValueError('Empty numeric selection')
        return (min if op=='min' else max)(float(r[spec['column']]) for r in rows)
    if op!='value' or len(rows)!=1:raise ValueError('Numeric selection must match exactly one row')
    value=rows[0][spec['column']]
    if isinstance(spec['expected'],bool):
        if value not in ['True','False','true','false']:raise ValueError('Nonboolean result')
        return value.lower()=='true'
    if isinstance(spec['expected'],(int,float)):return float(value)
    return value


def equal_value(actual,expected):
    if isinstance(expected,bool):return type(actual) is bool and actual==expected
    if isinstance(expected,(int,float)):
        return not isinstance(actual,bool) and isinstance(actual,(int,float)) and math.isfinite(actual) and math.isclose(actual,expected,rel_tol=1e-10,abs_tol=1e-12)
    return actual==expected
