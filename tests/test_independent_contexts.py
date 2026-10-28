import sys
from pathlib import Path
import unittest
import pandas as pd
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_independent_contexts import validate_counts,gene_rows
from audit_independent_metadata import parse


class IndependentContextGuards(unittest.TestCase):
    def test_normalized_fractional_values_are_not_rounded(self):
        with self.assertRaises(ValueError):validate_counts(pd.DataFrame({'a':[1.,2.25],'b':[2.,3.]},index=['x','y']))

    def test_invalid_counts_or_ambiguous_rows_fail(self):
        for x in [pd.DataFrame({'a':[1,2]},index=['x','x']),pd.DataFrame({'a':[1,-1]},index=['x','y']),pd.DataFrame({'a':[1,np.nan]},index=['x','y']),pd.DataFrame({'a':[0,0]},index=['x','y'])]:
            with self.assertRaises(ValueError):validate_counts(x)

    def test_gene_repeat_separation_uses_all_annotations(self):
        d=pd.DataFrame({'repName':[None,'L1',None],'repClass':[None,'LINE','LINE'],'repFamily':[None,'L1',None]})
        self.assertEqual(gene_rows(d).tolist(),[True,False,False])

    def test_valid_source_integer_counts_preserved(self):
        d=pd.DataFrame({'a':[0,2],'b':[3,0]},index=['x','y']);original=d.copy()
        validate_counts(d);pd.testing.assert_frame_equal(d,original)

    def test_metadata_keeps_age_conflicts_and_drops_contact_fields(self):
        series,samples=parse('^SERIES = GSE1\n!Series_geo_accession = GSE1\n!Series_contact_email = private@example.org\n^SAMPLE = GSM1\n!Sample_title = P10\n!Sample_characteristics_ch1 = age: P12\n!Sample_contact_phone = private\n')
        self.assertNotIn('Series_contact_email',series)
        self.assertEqual(samples['GSM1']['Sample_characteristics_ch1'],['age: P12'])
        self.assertNotIn('Sample_contact_phone',samples['GSM1'])


if __name__=='__main__':unittest.main()
