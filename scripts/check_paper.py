#!/usr/bin/env python3
"""Traceable manuscript checks, never a scientific-validity certification."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
from urllib.parse import unquote
from paper_support import ROOT,read_yaml,safe_path,digest,frontmatter,result_value,equal_value

STATES={'planned','draft','verified','superseded'}
MARKER=re.compile(r'<!--\s*claim:\s*([A-Z][A-Z0-9-]+)\s*-->')
CJK=re.compile(r'[\u3400-\u9fff]')
INFRA=['AGENTS.md','.agents/skills/paper-sync/SKILL.md','scripts/check_paper.py','scripts/paper_support.py',
       'scripts/export_paper.sh','scripts/export_paper.py','scripts/check_sync_readiness.py','docs/PAPER_WORKFLOW.md',
       'envs/requirements-paper.txt','envs/requirements-paper-check.txt','envs/paper_python.sh','.github/workflows/paper-check.yml']


class Checker:
    def __init__(self,root,submission=False):self.root=root;self.submission=submission;self.issues=[]
    def add(self,level,code,message):self.issues.append({'level':level,'code':code,'message':message})
    def error(self,code,message):self.add('ERROR',code,message)
    def warning(self,code,message):self.add('ERROR' if self.submission else 'WARNING',code,message)
    def mapping(self,relative):return read_yaml(safe_path(self.root,relative))
    def indexed(self,doc,key,prefix):
        rows=doc.get(key)
        if not isinstance(rows,list):raise ValueError(key+' must be a list')
        result={}
        for row in rows:
            if not isinstance(row,dict) or not re.fullmatch(prefix+r'-[A-Z0-9-]+',str(row.get('id',''))):raise ValueError('Invalid '+key+' ID')
            if row['id'] in result:self.error('DUPLICATE_ID',row['id'])
            result[row['id']]=row
            if row.get('status') not in STATES:self.error('STATUS',row['id']+' has invalid status')
        return result

    def infrastructure(self):
        for p in INFRA:
            if not safe_path(self.root,p).is_file():self.error('MISSING_INFRA',p)
        skill=safe_path(self.root,'.agents/skills/paper-sync/SKILL.md')
        if skill.exists():
            meta,_=frontmatter(skill)
            if meta.get('name')!='paper-sync' or not meta.get('description'):self.error('SKILL_METADATA','Missing skill identity/description')

    def evidence(self,rows):
        valid=set();source_sets={}
        for eid,e in rows.items():
            before=len([i for i in self.issues if i['level']=='ERROR']);sources=set()
            if e['status']=='verified':
                for field in ['verification_scope','recorded_on','version','configuration','reproduce','limitations']:
                    if not e.get(field):self.error('EVIDENCE_FIELD',eid+' missing '+field)
                if not e.get('files') and not e.get('snapshot'):self.error('EVIDENCE_EMPTY',eid+' has no source')
                if e.get('verification_scope')=='biological_validation':self.error('UNSUPPORTED_VALIDATION',eid+' cannot use machine checks as biological validation')
            for item in e.get('files',[]):
                try:
                    p=safe_path(self.root,item['path']);sources.add(item['path'])
                    if not p.is_file():raise ValueError('Missing source '+item['path'])
                    if not re.fullmatch(r'[0-9a-f]{64}',str(item.get('sha256',''))) or digest(p)!=item['sha256']:
                        raise ValueError('Changed source '+item['path'])
                except (OSError,ValueError,KeyError) as exc:self.error('EVIDENCE_FILE',eid+': '+str(exc))
            if e.get('snapshot'):
                try:
                    from check_project import verify_manifest
                    m=safe_path(self.root,e['snapshot']+'/manifest.json')
                    if digest(m)!=e.get('snapshot_manifest_sha256'):raise ValueError('Snapshot identity changed')
                    verify_manifest(self.root,e['snapshot'])
                    manifest=json.loads(m.read_text())
                    sources.update(e['snapshot']+'/'+r['path'] for r in manifest['files'])
                except (OSError,ValueError,KeyError) as exc:self.error('EVIDENCE_SNAPSHOT',eid+': '+str(exc))
            after=len([i for i in self.issues if i['level']=='ERROR'])
            if e['status']=='verified' and after==before:valid.add(eid)
            source_sets[eid]=sources
        return valid,source_sets

    def bibliography(self):
        try:
            import bibtexparser
            text=safe_path(self.root,'paper/references.bib').read_text()
            db=bibtexparser.loads(text)
            ids=[r['ID'] for r in db.entries]
            if len(ids)!=len(set(ids)):self.error('BIB_DUPLICATE','Duplicate BibTeX keys')
            if len(re.findall(r'(?m)^@(?:article|inproceedings|book|misc|techreport)\s*[{(]',text,re.I))!=len(db.entries):
                self.error('BIB_PARSE','Bibliography entries could not all be parsed')
            records=self.mapping('research/bibliography.yaml').get('references',[])
            if len({r['key'] for r in records})!=len(records):self.error('BIB_REGISTRY_DUPLICATE','Duplicate reference registry key')
            registry={r['key']:r for r in records};confirmed=set()
            for item in db.entries:
                key=item['ID'];record=registry.get(key,{})
                for field in ['title','author','year']:
                    if not item.get(field):self.error('BIB_FIELD',key+' missing '+field)
                doi=item.get('doi','')
                if doi and not re.fullmatch(r'10\.\d{4,9}/\S+',doi):self.error('BIB_DOI',key+' malformed DOI')
                if record.get('status')!='verified':self.warning('BIB_UNVERIFIED',key+' is not verified');continue
                if record.get('publication_status') not in ['published','preprint']:self.error('BIB_PUBLICATION',key+' missing publication status')
                if doi!=record.get('doi') or not record.get('verified_on') or not record.get('source_url'):
                    self.error('BIB_PROVENANCE',key+' missing/mismatched verification record');continue
                if record.get('metadata_file'):
                    p=safe_path(self.root,record['metadata_file'])
                    if digest(p)!=record.get('metadata_sha256'):self.error('BIB_METADATA',key+' metadata changed');continue
                    data=json.loads(p.read_text()).get('message',{})
                    if data.get('DOI','').lower()!=doi.lower():self.error('BIB_METADATA',key+' metadata DOI differs');continue
                    normalize=lambda s:re.sub(r'[^a-z0-9]','',s.lower())
                    if normalize(item['title'])!=normalize(data['title'][0]):self.error('BIB_TITLE',key+' title differs from source metadata');continue
                    year=str(data.get('published',data.get('issued',{}))['date-parts'][0][0])
                    authors=' and '.join((a.get('family','')+', '+a.get('given','')).strip(', ') for a in data.get('author',[]))
                    if str(item['year'])!=year or normalize(item.get('journal',''))!=normalize(data['container-title'][0]) or normalize(item['author'])!=normalize(authors):
                        self.error('BIB_IDENTITY',key+' author/year/venue differs from registered source metadata');continue
                confirmed.add(key)
            return set(ids),confirmed
        except ImportError:
            self.error('BIB_DEPENDENCY','Install bibtexparser from envs/requirements-paper-check.txt');return set(),set()

    def links(self,path,body):
        for image,label,target in re.findall(r'(!?)\[([^\]]*)\]\(([^)]+)\)',body):
            target=target.strip().strip('<>')
            if re.match(r'https?://|mailto:',target):continue
            if re.match(r'[A-Za-z][A-Za-z0-9+.-]*:',target):self.error('LINK_SCHEME',str(path.name)+' uses unsupported link');continue
            relative=unquote(target.split('#')[0])
            if not relative:continue
            p=(path.parent/relative).resolve()
            if not p.is_relative_to(self.root.resolve()) or not p.is_file():self.error('LINK_TARGET',str(path.name)+' broken/escaping link: '+target)
            if image and not label.strip():self.warning('IMAGE_ALT',str(path.name)+' has empty image caption/alt')

    def full(self):
        required=['paper/metadata.yaml','paper/outline.md','paper/references.bib','paper/figures/README.md','paper/templates/README.md',
            'research/progress.md','research/decisions.md','research/evidence.yaml','research/claims.yaml','research/terminology.yaml',
            'research/bibliography.yaml','research/paper-todo.md','research/submission-checklist.md','research/semantic-review.yaml']
        missing=[p for p in required if not safe_path(self.root,p).is_file()]
        for p in missing:self.error('MISSING_PRIVATE_FILE',p)
        if missing:return
        for folder in ['paper','research']:
            for p in sorted((self.root/folder).rglob('*')):
                if p.suffix in ['.yaml','.yml']:
                    try:read_yaml(p)
                    except (OSError,ValueError,TypeError) as exc:self.error('YAML_STRUCTURE',str(p.relative_to(self.root))+': '+str(exc))
        meta=self.mapping('paper/metadata.yaml')
        if meta.get('lang')!='en-US' or not meta.get('title') or CJK.search(str(meta.get('title'))):self.error('PAPER_LANGUAGE','Formal title/language must be English')
        for field in ['abstract','keywords']:
            if CJK.search(str(meta.get(field,''))):self.error('PAPER_LANGUAGE','Non-English '+field)
        chapters=meta.get('chapters',[])
        if not chapters or len(chapters)!=len(set(chapters)):self.error('CHAPTER_ORDER','Missing or repeated ordered chapters')
        actual={str(p.relative_to(self.root)) for p in (self.root/'paper/chapters').glob('*.md')}
        if set(chapters)!=actual:self.error('CHAPTER_INVENTORY','Ordered chapter list differs from actual chapters')
        evidence=self.indexed(self.mapping('research/evidence.yaml'),'evidence','E')
        claims=self.indexed(self.mapping('research/claims.yaml'),'claims','C')
        valid,sources=self.evidence(evidence);bib,confirmed=self.bibliography()
        termdoc=self.mapping('research/terminology.yaml');terms=termdoc.get('terms',[])
        if len({t['id'] for t in terms})!=len(terms):self.error('TERM_DUPLICATE','Duplicate terminology ID')
        blocks={};alltext=''
        for relative in chapters:
            p=safe_path(self.root,relative);fm,body=frontmatter(p);alltext+='\n'+body
            if fm.get('status') not in STATES or not fm.get('purpose') or not fm.get('id'):self.error('CHAPTER_METADATA',relative+' missing purpose/status/ID')
            if fm.get('status')!='verified':self.warning('CHAPTER_DRAFT',relative+' remains '+str(fm.get('status')))
            if CJK.search(body):self.error('PAPER_LANGUAGE',relative+' contains CJK prose/caption')
            if len(body.strip())<100 or len(re.findall(r'^# ',body,re.M))!=1:self.error('CHAPTER_STRUCTURE',relative+' needs one main heading and substantive content')
            levels=[len(x) for x in re.findall(r'^(#{1,6}) ',body,re.M)]
            if any(b>a+1 for a,b in zip(levels,levels[1:])):self.error('HEADING_JUMP',relative)
            if re.search(r'\bTODO\b|\bTBD\b|\[VERIFY\]|待验证',body,re.I):self.warning('UNRESOLVED_MARKER',relative)
            self.links(p,body)
            for key in re.findall(r'(?<![\w@])@([A-Za-z][A-Za-z0-9_:-]*)',body):
                if key not in bib:self.error('CITATION_MISSING',relative+': '+key)
                elif key not in confirmed:self.warning('CITATION_UNVERIFIED',relative+': '+key)
            matches=list(MARKER.finditer(body))
            for i,m in enumerate(matches):
                cid=m[1]
                if cid not in claims:self.error('CLAIM_MARKER',relative+': unknown '+cid)
                blocks.setdefault(cid,[]).append((relative,body[m.end():matches[i+1].start() if i+1<len(matches) else len(body)]))
            for t in terms:
                for banned in t.get('prohibited_phrases',[]):
                    if re.search(r'(?<!\w)'+re.escape(banned)+r'(?!\w)',body,re.I):self.warning('TERM_REVIEW',relative+' uses '+banned)
        for cid,c in claims.items():
            ids=c.get('evidence_ids',[])
            if any(e not in evidence for e in ids):self.error('CLAIM_EVIDENCE',cid+' references missing evidence')
            if c['status']=='verified' and (not ids or not set(ids)<=valid):self.error('CONFIRMED_UNSUPPORTED',cid+' lacks current scoped verified evidence')
            if c['status']=='verified':
                for _,block in blocks.get(cid,[]):
                    if any(key not in confirmed for key in re.findall(r'(?<![\w@])@([A-Za-z][A-Za-z0-9_:-]*)',block)):
                        self.error('CONFIRMED_CITATION',cid+' uses an unverified citation')
            if c.get('kind') in ['interpretation','hypothesis'] and c['status']=='verified':self.error('INFERENCE_STATUS',cid+' interpretation cannot be promoted by structural checks')
            if c.get('requires_human_review'):self.warning('HUMAN_REVIEW',cid+' requires scientific review')
            if c['status'] in ['planned','superseded'] and blocks.get(cid):self.error('INACTIVE_CLAIM',cid+' used in prose')
            locations={p for p,b in blocks.get(cid,[])}
            if c['status'] not in ['planned','superseded'] and set(c.get('chapters',[]))!=locations:self.error('CLAIM_LOCATION',cid+' registry/prose locations differ')
            allowed=set().union(*(sources.get(e,set()) for e in ids)) if ids else set()
            for a in c.get('assertions',[]):
                try:
                    if a['source'] not in allowed:raise ValueError('Numeric source is not linked evidence')
                    actual_value=result_value(self.root,a)
                    if not equal_value(actual_value,a['expected']):raise ValueError('Structured result differs from expected value')
                    matching=[b for p,b in blocks.get(cid,[]) if p==a['chapter']]
                    if not matching or not any(a['text'] in b for b in matching):raise ValueError('Expected numeric text absent from claim block')
                except (OSError,ValueError,KeyError,TypeError) as exc:self.error('NUMERIC_ASSERTION',cid+': '+str(exc))
        for relative in ['paper/outline.md','paper/figures/README.md','paper/templates/README.md','research/progress.md','research/decisions.md','research/paper-todo.md','research/submission-checklist.md']:
            p=safe_path(self.root,relative);self.links(p,p.read_text())
        review=self.mapping('research/semantic-review.yaml')
        if review.get('reviewer')!='Codex' or review.get('scientific_validation') is not False:self.error('SEMANTIC_SCOPE','Semantic review must not claim biological validation')
        for relative,checksum in review.get('reviewed_files',{}).items():
            if digest(safe_path(self.root,relative))!=checksum:self.warning('SEMANTIC_STALE','Review stale for '+relative)
        if not set(chapters)<=set(review.get('reviewed_files',{})):self.warning('SEMANTIC_COVERAGE','Not all chapters have a hash-bound semantic review')
        pending=sum(p.read_text().count('- [ ]') for p in [self.root/'research/paper-todo.md',self.root/'research/submission-checklist.md'])
        if pending:self.warning('OPEN_ITEMS',str(pending)+' tracked manuscript/submission items remain')
        if not meta.get('authors') or not meta.get('target_venue'):self.warning('SUBMISSION_METADATA','Authorship and target venue remain unconfirmed')
        self.add('INFO','CHECK_BOUNDARY','Checks establish structure, traceability and selected numeric agreement; not scientific truth, novelty or submission quality.')
        self.add('INFO','COVERAGE',f'{len(chapters)} chapters, {len(evidence)} evidence records, {len(claims)} claims, {len(bib)} references checked')


def check(root,mode='full',submission=False):
    c=Checker(Path(root).resolve(),submission)
    try:
        c.infrastructure()
        if mode=='full':c.full()
        else:c.add('INFO','INFRA_SCOPE','Infrastructure mode does not validate a private manuscript, evidence or submission readiness.')
    except (OSError,ValueError,KeyError,TypeError) as exc:c.error('INVALID_DOCUMENT',str(exc))
    return c.issues


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--infrastructure-only',action='store_true');p.add_argument('--submission-ready',action='store_true');p.add_argument('--report',type=Path);args=p.parse_args()
    if args.infrastructure_only and args.submission_ready:p.error('Infrastructure mode cannot certify submission readiness')
    issues=check(args.root,'infrastructure' if args.infrastructure_only else 'full',args.submission_ready)
    for row in issues:print(f"{row['level']} {row['code']}: {row['message']}")
    counts=Counter(r['level'] for r in issues);print('Summary:',json.dumps(dict(counts),sort_keys=True))
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps({'mode':'infrastructure' if args.infrastructure_only else 'full','counts':dict(counts),'issues':issues},indent=2)+'\n')
    return 1 if counts['ERROR'] else 0


if __name__=='__main__':raise SystemExit(main())
