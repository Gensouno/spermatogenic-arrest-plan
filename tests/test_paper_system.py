"""Synthetic regression fixtures only; no unpublished research is embedded here."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from paper_support import yaml,read_yaml,digest,result_value,safe_path
from check_paper import check,INFRA
from check_sync_readiness import private_path,inspect
from export_paper import rewrite_docx,W,q


def fixture(root):
    def text(path,content):
        p=root/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(content)
    def y(path,value):text(path,yaml.safe_dump(value,sort_keys=False))
    for name in INFRA:text(name,'Infrastructure fixture\n')
    text('.agents/skills/paper-sync/SKILL.md','---\nname: paper-sync\ndescription: Synthetic test skill\n---\nFixture only.\n')
    chapter='paper/chapters/results.md'
    text(chapter,'---\nid: results\nstatus: draft\npurpose: Test evidence checking\n---\n# Results\n\n<!-- claim: C-TEST -->\nThis synthetic regression fixture reports 3 observations and does not represent a real research result. It exists only to test traceability checks.\n')
    text('result.csv','group,value\na,3\n')
    y('paper/metadata.yaml',{'lang':'en-US','title':'Synthetic checker fixture','chapters':[chapter],'authors':[],'target_venue':None})
    for name in ['paper/outline.md','paper/figures/README.md','paper/templates/README.md','research/progress.md',
                 'research/decisions.md','research/paper-todo.md','research/submission-checklist.md']:text(name,'# Synthetic fixture\n')
    y('research/evidence.yaml',{'evidence':[{'id':'E-TEST','status':'verified','verification_scope':'synthetic_numeric_test',
        'recorded_on':'2026-10-27','version':'fixture','configuration':'fixture','reproduce':['synthetic fixture'],
        'limitations':['Not a scientific result'],'files':[{'path':'result.csv','sha256':digest(root/'result.csv')}]}]})
    y('research/claims.yaml',{'claims':[{'id':'C-TEST','kind':'computational_observation','status':'verified','evidence_ids':['E-TEST'],
        'chapters':[chapter],'requires_human_review':False,'assertions':[{'source':'result.csv','format':'csv','where':{'group':'a'},
        'column':'value','expected':3,'chapter':chapter,'text':'3 observations'}]}]})
    y('research/terminology.yaml',{'terms':[]})
    text('paper/references.bib','@article{Law2014, title={voom: precision weights unlock linear model analysis tools for RNA-seq read counts}, author={Law, Charity W. and Chen, Yunshun and Shi, Wei and Smyth, Gordon K.}, year={2014}, journal={Genome Biology}, doi={10.1186/gb-2014-15-2-r29}}\n')
    y('research/bibliography.yaml',{'references':[{'key':'Law2014','status':'verified','publication_status':'published',
        'doi':'10.1186/gb-2014-15-2-r29','verified_on':'2026-10-27','source_url':'https://doi.org/10.1186/gb-2014-15-2-r29'}]})
    y('research/semantic-review.yaml',{'reviewer':'Codex','scientific_validation':False,'reviewed_files':{chapter:digest(root/chapter)}})
    return chapter


class PaperChecksTests(unittest.TestCase):
    def setUp(self):self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name);self.chapter=fixture(self.root)
    def tearDown(self):self.temp.cleanup()
    def errors(self):return {r['code'] for r in check(self.root) if r['level']=='ERROR'}
    def modify_yaml(self,name,change):
        p=self.root/name;d=read_yaml(p);change(d);p.write_text(yaml.safe_dump(d,sort_keys=False))

    def test_valid_draft_passes_but_is_not_submission_ready(self):
        self.assertEqual(self.errors(),set())
        self.assertTrue(any(r['level']=='ERROR' for r in check(self.root,submission=True)))

    def test_duplicate_yaml_key_rejected(self):
        p=self.root/'duplicate.yaml';p.write_text('a: 1\na: 2\n')
        with self.assertRaises(ValueError):read_yaml(p)

    def test_malformed_optional_yaml_produces_structured_error(self):
        (self.root/'research/other.yaml').write_text('unfinished: [\n')
        self.assertIn('YAML_STRUCTURE',self.errors())

    def test_changed_result_invalidates_evidence_and_confirmed_claim(self):
        (self.root/'result.csv').write_text('group,value\na,4\n')
        self.assertTrue({'EVIDENCE_FILE','CONFIRMED_UNSUPPORTED','NUMERIC_ASSERTION'}<=self.errors())

    def test_ambiguous_numeric_selection_is_not_silently_first_row(self):
        (self.root/'result.csv').write_text('group,value\na,3\na,3\n')
        with self.assertRaises(ValueError):result_value(self.root,{'source':'result.csv','format':'csv','where':{'group':'a'},'column':'value','expected':3})

    def test_missing_evidence_and_duplicate_claim_ids_fail(self):
        self.modify_yaml('research/claims.yaml',lambda d:d['claims'][0].update(evidence_ids=['E-MISSING']))
        self.assertIn('CLAIM_EVIDENCE',self.errors())
        self.modify_yaml('research/claims.yaml',lambda d:d['claims'].append(d['claims'][0].copy()))
        self.assertIn('DUPLICATE_ID',self.errors())

    def test_draft_evidence_cannot_support_confirmed_claim(self):
        self.modify_yaml('research/evidence.yaml',lambda d:d['evidence'][0].update(status='draft'))
        self.assertIn('CONFIRMED_UNSUPPORTED',self.errors())

    def test_unverified_citation_cannot_support_confirmed_claim(self):
        p=self.root/self.chapter;p.write_text(p.read_text()+'\nSupported by [@Law2014].\n')
        self.modify_yaml('research/bibliography.yaml',lambda d:d['references'][0].update(status='draft'))
        self.assertIn('CONFIRMED_CITATION',self.errors())

    def test_numeric_prose_drift_detected(self):
        p=self.root/self.chapter;p.write_text(p.read_text().replace('3 observations','4 observations'))
        self.assertIn('NUMERIC_ASSERTION',self.errors())

    def test_english_caption_and_missing_reference_checked(self):
        p=self.root/self.chapter;p.write_text(p.read_text()+'\n中文图题 [@MissingSource]\n')
        self.assertTrue({'PAPER_LANGUAGE','CITATION_MISSING'}<=self.errors())

    def test_link_and_source_path_escape_rejected(self):
        p=self.root/self.chapter;p.write_text(p.read_text()+'\n![Outside](../../../outside.png)\n')
        self.assertIn('LINK_TARGET',self.errors())
        with self.assertRaises(ValueError):safe_path(self.root,'../outside')

    def test_infrastructure_mode_does_not_require_or_claim_private_paper_check(self):
        (self.root/'paper/metadata.yaml').unlink()
        result=check(self.root,mode='infrastructure')
        self.assertFalse(any(r['level']=='ERROR' for r in result))
        self.assertIn('INFRA_SCOPE',{r['code'] for r in result})
        self.assertIn('MISSING_PRIVATE_FILE',self.errors())

    def test_private_and_secret_material_block_release_without_echoing_secret(self):
        self.assertTrue(private_path('paper/chapters/results.md'))
        self.assertTrue(private_path('research/claims.yaml'))
        self.assertFalse(private_path('scripts/check_paper.py'))
        secret=('gh'+'p_'+'A'*40).encode();(self.root/'unsafe.py').write_bytes(secret)
        def fakegit(root,*args):
            if args==('ls-files','--cached','--others','--exclude-standard','-z'):return b'paper/chapters/results.md\0unsafe.py\0'
            if args==('ls-files','-z'):return b'paper/chapters/results.md\0'
            if args==('status','--porcelain=v1','-z'):return b'?? unsafe.py\0'
            if args==('rev-parse','HEAD'):return b'a'*40+b'\n'
            if args==('branch','--show-current'):return b'main\n'
            if args==('remote',):return b''
            raise AssertionError(args)
        with patch('check_sync_readiness.git',side_effect=fakegit):report=inspect(self.root)
        self.assertEqual({r['code'] for r in report['issues']},{'PROTECTED_PATH','POSSIBLE_SECRET'})
        self.assertNotIn(secret.decode(),json.dumps(report))

    def test_word_reference_preserves_pandoc_style_prefix_and_page_size(self):
        source=self.root/'source.docx';target=self.root/'reference.docx'
        with ZipFile(source,'w') as z:
            z.writestr('word/styles.xml',f'<w:styles xmlns:w="{W}"><w:style w:type="paragraph" w:styleId="Normal"/></w:styles>')
            z.writestr('word/document.xml',f'<w:document xmlns:w="{W}"><w:body><w:sectPr/></w:body></w:document>')
        rewrite_docx(source,target,reference=True)
        with ZipFile(target) as z:
            styles=z.read('word/styles.xml').decode()
            # This is how Pandoc extends a reference file; the injected prefix
            # must remain bound even after our style edits.
            extended=styles.replace('</w:styles>','<w:style w:type="paragraph" w:styleId="Injected"/></w:styles>')
            self.assertEqual(len(ET.fromstring(extended).findall(q('style'))),2)
            page=ET.fromstring(z.read('word/document.xml')).find('.//'+q('pgSz'))
            self.assertEqual((page.get(q('w')),page.get(q('h'))),('11906','16838'))
        # A supplied reference's page size must survive final formatting.
        custom=self.root/'custom.docx';preserved=self.root/'preserved.docx'
        with ZipFile(custom,'w') as z:
            z.writestr('word/document.xml',f'<w:document xmlns:w="{W}"><w:body><w:sectPr><w:pgSz w:w="12240" w:h="15840"/></w:sectPr></w:body></w:document>')
        rewrite_docx(custom,preserved,table_style=False)
        with ZipFile(preserved) as z:
            page=ET.fromstring(z.read('word/document.xml')).find('.//'+q('pgSz'))
            self.assertEqual((page.get(q('w')),page.get(q('h'))),('12240','15840'))


if __name__=='__main__':unittest.main()
