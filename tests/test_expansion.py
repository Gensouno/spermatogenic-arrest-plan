import sys
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
sys.path.insert(0, str(ROOT / 'scripts'))
from prepare_expansion import validate_counts, featurecounts, parse_single_table
from discover_expansion import parse_soft


class ExpansionTests(unittest.TestCase):
    def test_reject_fractional_abundance_without_rounding(self):
        with self.assertRaises(ValueError):
            validate_counts(pd.DataFrame({'s1': [1, 1.5]}))

    def test_reject_duplicate_gene_identifiers(self):
        with self.assertRaises(ValueError):
            validate_counts(pd.DataFrame({'s1': [1, 2]}, index=['a', 'a']))

    def test_reject_nonfinite_or_empty_library(self):
        for x in [pd.DataFrame({'a': [1, np.nan]}), pd.DataFrame({'a': [0, 0]})]:
            with self.assertRaises(ValueError):
                validate_counts(x)

    def test_featurecounts_does_not_import_length_as_counts(self):
        value = featurecounts(b'# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tsample.bam\na\t1\t1\t101\t+\t100\t4\n')
        self.assertEqual(value['a'], 4)
        self.assertEqual(value.name, 'sample.bam')

    def test_comma_and_tab_sample_values_are_identical(self):
        a = parse_single_table(b'gene,sample\na,3\nb,4\n')
        b = parse_single_table(b'gene\tsample\na\t3\nb\t4\n')
        pd.testing.assert_series_equal(a, b)

    def test_soft_preserves_genotype_and_stops_at_platform(self):
        text = '^SERIES = GSE1\n!Series_title = source\n^SAMPLE = GSM1\n!Sample_characteristics_ch1 = genotype: Pramef12+/-\n^PLATFORM = GPL1\n!Sample_title = must_not_leak\n'
        series, samples = parse_soft(text)
        self.assertEqual(samples['GSM1']['Sample_characteristics_ch1'], ['genotype: Pramef12+/-'])
        self.assertNotIn('Sample_title', samples['GSM1'])


if __name__ == '__main__':
    unittest.main()
