import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'scripts')]
from prepare_age_reference import barcode_assignments,parse_count_row,read_author_table
from diagnose_reference_transfer import gene_folds,train_fit
from check_project import safe_path,verify_manifest


class ReferenceTransferTests(unittest.TestCase):
    def test_strict_ooxml_sample_table_without_rewriting_source(self):
        namespace='http://purl.oclc.org/ooxml/spreadsheetml/main'
        strings=['sample','Conditions']+[f'library{i}' for i in range(8)]+['period']
        shared=f'<sst xmlns="{namespace}">'+''.join('<si><t>'+x+'</t></si>' for x in strings)+'</sst>'
        rows=['<row><c t="s"><v>0</v></c><c t="s"><v>1</v></c></row>']
        rows += [f'<row><c t="s"><v>{i+2}</v></c><c t="s"><v>10</v></c></row>' for i in range(8)]
        sheet=f'<worksheet xmlns="{namespace}"><sheetData>'+''.join(rows)+'</sheetData></worksheet>'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source.xlsx'
            with zipfile.ZipFile(p,'w') as z:
                z.writestr('xl/sharedStrings.xml',shared);z.writestr('xl/worksheets/sheet1.xml',sheet)
            before=p.read_bytes();frame=read_author_table(p)
            self.assertEqual(frame.shape,(8,2));self.assertEqual(frame.iloc[0,0],'library0')
            self.assertEqual(p.read_bytes(),before)

    def test_count_row_rejects_fractional_and_truncated_fields(self):
        for line in ['"A",1,2.5','"A",1,2,garbage','"A",1']:
            with self.assertRaises(ValueError):parse_count_row(line,2)
        gene,x=parse_count_row('"A",1,2\n',2)
        self.assertEqual(gene,'A');np.testing.assert_array_equal(x,[1,2])

    def test_barcode_source_must_be_known_and_complete(self):
        with self.assertRaises(ValueError):barcode_assignments(['unknown_AAAA'])
        with self.assertRaises(ValueError):barcode_assignments(['PND06_AAAA','PND06_AAAA'])

    def test_folds_are_deterministic_and_never_overlap(self):
        a=gene_folds([str(i) for i in range(25)])
        self.assertEqual(a,gene_folds([str(i) for i in range(25)]))
        self.assertEqual(set(a.values()),set(range(5)))
        with self.assertRaises(ValueError):gene_folds(['A']*5)

    def test_heldout_target_values_do_not_change_training_coefficients(self):
        r=np.array([[1.,0],[0,1],[1,1],[2,1]])
        y=r@np.array([.3,.7]);train=np.array([True,True,False,False])
        before=train_fit(r,y,np.ones(4),train)
        y[~train]=[999,999]
        np.testing.assert_allclose(before,train_fit(r,y,np.ones(4),train))

    def test_manifest_catches_content_and_source_changes(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);base=root/'snapshot';base.mkdir()
            (base/'x.csv').write_text('a\n1\n');(root/'source.py').write_text('pass\n')
            h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
            m={'files':[{'path':'x.csv','bytes':4,'sha256':h(base/'x.csv')}],
               'source_sha256':{'source.py':h(root/'source.py')}}
            (base/'manifest.json').write_text(json.dumps(m))
            self.assertEqual(verify_manifest(root,'snapshot')['files'],1)
            (base/'nested').mkdir();(base/'nested/manifest.json').write_text('{}')
            with self.assertRaises(ValueError):verify_manifest(root,'snapshot')
            (base/'nested/manifest.json').unlink()
            (root/'source.py').write_text('changed\n')
            with self.assertRaises(ValueError):verify_manifest(root,'snapshot')

    def test_manifest_paths_cannot_escape(self):
        for name in ['../outside','/tmp/outside']:
            with self.assertRaises(ValueError):safe_path(ROOT,name)


if __name__=='__main__':unittest.main()
