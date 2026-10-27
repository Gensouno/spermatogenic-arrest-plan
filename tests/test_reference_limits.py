"""Reference identity, exon accounting, source aliases and holdout failure modes."""
import gzip
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'code'),str(ROOT/'scripts')]
from diagnose_reference_limits import merged_length,exon_lengths,extended_folds
from diagnose_reference_transfer import gene_folds
from check_chen_transfer import candidate_weights,specificity_panel
from check_chen_library_holdout import leave_library_profile
from prepare_chen_reference import read_umi,broad_stage
import fetch_chen_umi


class ReferenceLimitsTests(unittest.TestCase):
    def test_inclusive_exon_union_not_transcript_sum(self):
        self.assertEqual(merged_length([(1,10),(4,15),(20,20),(21,23),(1,10)]),19)
        self.assertEqual(merged_length([]),0)
        for intervals in [[(0,2)],[(5,4)]]:
            with self.assertRaises(ValueError):merged_length(intervals)

    def test_length_parser_rejects_multiple_chromosome_gene(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'gtf.gz'
            with gzip.open(path,'wt') as f:f.write('1\tx\texon\t1\t10\t.\t+\t.\tgene_id "ENSMUSG1";\n2\tx\texon\t1\t10\t.\t+\t.\tgene_id "ENSMUSG1";\n')
            with self.assertRaises(ValueError):exon_lengths(path,pd.Series({'ENSMUSG1':'Gene'}))

    def test_old_folds_retained_when_new_genes_added(self):
        old=list('ABCDEFGHIJ');new=old+['AA','K','ZZ']
        result=extended_folds(old,new)
        self.assertEqual({g:result[g] for g in old},gene_folds(sorted(old)))
        self.assertEqual(result,extended_folds(old,list(reversed(new))))

    def test_extra_components_retain_broad_mass(self):
        columns=['Somatic','GC1',*[f'GC{i}' for i in range(2,13)],'Chen_L','Chen_Z','Chen_BS']
        p,b=candidate_weights(np.ones(len(columns)),columns)
        self.assertAlmostEqual(b.sum(),1)
        self.assertAlmostEqual(b['Spermatogonia'],2/len(columns))
        self.assertAlmostEqual(b['Meiotic'],9/len(columns))
        with self.assertRaises(ValueError):candidate_weights(np.ones(1),['unknown'])

    def test_source_marker_panel_invariant_to_row_order(self):
        chen=pd.DataFrame({'L':[5,0,1,1],'BS':[0,5,1,1]},index=['Lgene','BSgene','TieA','TieB'])
        green=pd.DataFrame({'Somatic':[0,0,4,4],'GC1':[5,5,1,1]},index=chen.index)
        first,_=specificity_panel(chen,green,per_component=1)
        second,_=specificity_panel(chen.iloc[::-1],green.iloc[::-1],per_component=1)
        self.assertEqual(first,second);self.assertEqual(first,{'Lgene','BSgene','TieA'})

    def test_library_mean_recombined_by_cells(self):
        e=pd.DataFrame({'profile':['a','b','c'],'author_type':['L','L','L'],'library':['A','B','C'],
            'weighting':['equal_cell_mean']*3,'cells':[1,5,4]})
        x=pd.DataFrame({'a':[100.],'b':[200.],'c':[900.]},index=['Gene'])
        v,n,libs=leave_library_profile(x,e,'L','C')
        self.assertAlmostEqual(v.iloc[0],1100/6);self.assertEqual(n,6);self.assertEqual(libs,['A','B'])
        with self.assertRaises(ValueError):leave_library_profile(x,e,'Z','C')

    def test_explicit_cell_alias_preserves_author_ID(self):
        j=pd.DataFrame({'cell':['SpSC24_In_sc.3'],'raw_file':[np.nan],'processed_file':[np.nan],'barcode':[np.nan]})
        b=pd.DataFrame({'cell':['SpSC24_ln_sc.3'],'raw_file':['SpSC24'],'processed_file':['SpSC24_UMIcount.txt'],'barcode':['ACGT']})
        r=fetch_chen_umi.reconcile_id_aliases(j,b)
        self.assertEqual(r.iloc[0].cell,'SpSC24_In_sc.3');self.assertEqual(r.iloc[0].matrix_cell,'SpSC24_ln_sc.3')
        self.assertTrue(r.iloc[0].explicit_ID_alias_used)
        with self.assertRaises(ValueError):fetch_chen_umi.reconcile_id_aliases(j,pd.concat([b,b]))

    def test_alias_not_generalized_to_unlisted_library(self):
        j=pd.DataFrame({'cell':['SpSC99_In_sc.3'],'raw_file':[np.nan],'processed_file':[np.nan],'barcode':[np.nan]})
        b=pd.DataFrame({'cell':['SpSC99_ln_sc.3'],'raw_file':['SpSC99'],'processed_file':['x'],'barcode':['ACGT']})
        r=fetch_chen_umi.reconcile_id_aliases(j,b)
        self.assertFalse(r.iloc[0].explicit_ID_alias_used);self.assertEqual(r.iloc[0].matrix_cell,'SpSC99_In_sc.3')

    def test_fractional_umi_and_duplicate_gene_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'umi.tsv'
            for text in ['Gene\tc\nA\t1.5\n','Gene\tc\nA\t1\nA\t2\n']:
                path.write_text(text)
                with self.assertRaises(ValueError):read_umi(path)

    def test_preleptotene_not_forced_to_meiotic(self):
        for stage in ['G1','ePL','mPL','lPL']:self.assertIsNone(broad_stage(stage))
        self.assertEqual(broad_stage('L'),'Meiotic');self.assertEqual(broad_stage('BS'),'Spermatogonia')
        with self.assertRaises(ValueError):broad_stage('NewCell')

    def test_derived_workbook_restored_offline_from_bundle(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(fetch_chen_umi,'ROOT',Path(temp)):
            bundle=Path(temp)/'bundle.zip'
            with zipfile.ZipFile(bundle,'w') as z:z.writestr('table',b'abc')
            row={'path':'table','source_bundle':'bundle.zip','member':'table','bytes':3,'sha256':hashlib.sha256(b'abc').hexdigest()}
            fetch_chen_umi.restore_asset(row,False)
            self.assertEqual((Path(temp)/'table').read_bytes(),b'abc')


if __name__=='__main__':unittest.main()
