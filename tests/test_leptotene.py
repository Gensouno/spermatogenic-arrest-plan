"""Meaningful source-count and identifier failure modes for stage diagnostics."""
import gzip
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'scripts')]
from diagnose_leptotene import checked_counts,identifier_map,read_htseq,split_htseq_rows,paper_evidence
import fetch_leptotene


class LeptoteneTests(unittest.TestCase):
    def test_summary_counts_never_enter_gene_coverage(self):
        t=pd.DataFrame({'sample':[10,20,900,70]},index=['Gene1','Gene2','__not_aligned','__no_feature'])
        genes,summary=split_htseq_rows(t)
        self.assertEqual(genes['sample'].sum(),30)
        self.assertEqual(summary['sample'].sum(),970)
        self.assertEqual(genes.loc[['Gene1'],'sample'].sum()/genes['sample'].sum(),1/3)

    def test_unknown_summary_not_silently_discarded(self):
        with self.assertRaises(ValueError):split_htseq_rows(pd.DataFrame({'a':[1,2]},index=['Gene1','__unknown']))

    def test_bad_counts_refused(self):
        for values in [[1,-1],[1,.5],[1,np.nan],[0,0]]:
            with self.assertRaises(ValueError):checked_counts(pd.DataFrame({'a':values},index=['Gene1','Gene2']))
        with self.assertRaises(ValueError):checked_counts(pd.DataFrame({'a':[1,2]},index=['Gene','Gene']))

    def test_mapping_discards_all_ambiguous_symbols(self):
        genes=pd.DataFrame({'ID':['ENSMUSG1','ENSMUSG2','ENSMUSG3','ENSG4'],'Symbol':['dup','dup','ok','human']})
        mapping=identifier_map(genes)
        self.assertEqual(mapping.to_dict(),{'ENSMUSG3':'ok'})

    def test_htseq_gene_mapping_and_summary_audit(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'sample.gz'
            with gzip.open(path,'wt') as out:out.write('ENSMUSG1\t10\nENSMUSG2\t30\n__no_feature\t999\n')
            vector,audit=read_htseq(path,pd.Series({'ENSMUSG1':'Gene1'}))
            self.assertEqual(vector.to_dict(),{'Gene1':10})
            self.assertEqual(audit['mapping_expression_mass'],.25)
            self.assertEqual(audit['gene_count_total'],40)
            self.assertEqual(audit['excluded_HTSeq_summary_counts'],{'__no_feature':999})

    def test_versioned_ids_not_silently_converted(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'sample.gz'
            with gzip.open(path,'wt') as out:out.write('ENSMUSG1.2\t10\n')
            with self.assertRaises(ValueError):read_htseq(path,pd.Series({'ENSMUSG1':'Gene1'}))

    def test_wrong_primary_paper_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'paper.xml';path.write_text('<article><article-id pub-id-type="doi">10.1/wrong</article-id></article>')
            with self.assertRaises(ValueError):paper_evidence(path)

    def test_frozen_missing_input_stops_offline(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(fetch_leptotene,'ROOT',Path(temp)):
            with self.assertRaises(FileNotFoundError):fetch_leptotene.verify_asset({'path':'source'},False)

    def test_restored_input_checked_before_promotion(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(fetch_leptotene,'ROOT',Path(temp)):
            row={'path':'source','bytes':3,'url':'https://example.invalid/source','sha256':hashlib.sha256(b'yes').hexdigest()}
            def bad_fetch(path,*args):path.write_bytes(b'bad')
            with patch.object(fetch_leptotene,'fetch',side_effect=bad_fetch),self.assertRaises(ValueError):fetch_leptotene.verify_asset(row,True)
            self.assertFalse((Path(temp)/'source').exists())
            def good_fetch(path,*args):path.write_bytes(b'yes')
            with patch.object(fetch_leptotene,'fetch',side_effect=good_fetch):fetch_leptotene.verify_asset(row,True)
            self.assertEqual((Path(temp)/'source').read_bytes(),b'yes')

    def test_different_frozen_input_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(fetch_leptotene,'ROOT',Path(temp)):
            p=Path(temp)/'source';p.write_bytes(b'bad')
            row={'path':'source','bytes':3,'sha256':hashlib.sha256(b'yes').hexdigest()}
            with patch.object(fetch_leptotene,'fetch',side_effect=AssertionError('Do not download')),self.assertRaises(ValueError):fetch_leptotene.verify_asset(row,True)
            self.assertEqual(p.read_bytes(),b'bad')


if __name__=='__main__':unittest.main()
