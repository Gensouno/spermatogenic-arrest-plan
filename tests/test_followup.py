import sys
from pathlib import Path
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from prepare_followup import lrh1_samples, parse_fpkm
from summarize_followup import unique_symbol_effects


def source_fixture():
    rows = {}
    for i in range(16):
        age = 'P1' if i<8 else 'P3';control = i%8<4
        driver = 'Nanos3Cre' if i<8 else 'DhhCre'
        genotype = 'Lrh1 flox/-' if i<8 else 'Lrh1 flox/flox'
        if not control:
            genotype += '; '+('Nanos3-Cre' if i<8 else 'Dhh-Cre')
        rows[f'GSM{5455535+i}'] = {'Sample_title': [f'test [{age.lower()}_testis_{"noCre" if control else driver}_rep{i%4+1}]'],
            'Sample_characteristics_ch1': ['tissue: testis', 'age: Postnatal Day '+age[1:], 'genotype: '+genotype],
            'Sample_organism_ch1': ['Mus musculus']}
    return rows


class FollowupTests(unittest.TestCase):
    def test_alias_does_not_rename_floxed_controls_WT(self):
        frame = lrh1_samples(source_fixture())
        self.assertEqual(len(frame), 16)
        self.assertEqual(frame.cohort.nunique(), 2)
        self.assertEqual(frame.iloc[4].author_sample_name, 'p1_testis_NanosCre_rep1')
        self.assertEqual(frame.iloc[0].genotype, 'Lrh1 flox/-')

    def test_source_genotype_mismatch_stops_import(self):
        source = source_fixture()
        source['GSM5455535']['Sample_characteristics_ch1'][-1] = 'genotype: WT'
        with self.assertRaises(ValueError): lrh1_samples(source)

    def test_utf16_FPKM_is_explicitly_identified_not_counts(self):
        v = parse_fpkm('gene_id\tZ-01.FPKM\nENSMUST000001\t1.5\n'.encode('utf-16'))
        self.assertEqual(v.iloc[0], 1.5)
        with self.assertRaises(ValueError):
            parse_fpkm('gene_id\tcounts\nENSMUST000001\t2\n'.encode('utf-16'))

    def test_ambiguous_symbols_are_all_excluded_without_effect_selection(self):
        frame = pd.DataFrame({'symbol': ['A','A','B',None], 'logFC': [10,-1,2,99]})
        self.assertEqual(unique_symbol_effects(frame).to_dict(), {'B': 2})


if __name__ == '__main__':
    unittest.main()
