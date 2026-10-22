"""Fast synthetic checks; real source-table checks run when available locally."""
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
from check_stage_reference import collapse_centroids, fit, read_centroids, FILENAME


class StageReferenceTests(unittest.TestCase):
    def test_collapse_uses_cell_numbers_not_equal_cluster_weights(self):
        fine = pd.DataFrame({f'GC{i}': [float(i), 100 - i] for i in range(1, 13)}, index=['a', 'b'])
        cells = pd.Series({f'GC{i}': i for i in range(1, 13)})
        broad = collapse_centroids(fine, cells)
        self.assertAlmostEqual(broad.loc['a', 'Meiotic'], sum(i*i for i in range(2, 9)) / sum(range(2, 9)))
        np.testing.assert_allclose(broad.sum(), 100)
        np.testing.assert_allclose(broad.Spermatogonia, fine.GC1)

    def test_exact_known_linear_mixture(self):
        s = np.random.default_rng(1).gamma(2, 10, (100, 4))
        truth = np.array([.1, .2, .3, .4])
        np.testing.assert_allclose(fit(s, s @ truth, np.ones(100)), truth, atol=1e-12)

    def test_gene_holdout_does_not_use_test_target(self):
        s = np.random.default_rng(2).gamma(2, 10, (100, 4))
        y = s @ np.array([.1, .2, .3, .4])
        train = np.arange(100) % 5 != 0
        before = fit(s[train], y[train], np.ones(train.sum()))
        y[~train] *= 1000
        after = fit(s[train], y[train], np.ones(train.sum()))
        np.testing.assert_array_equal(before, after)

    def test_empty_target_rejected(self):
        with self.assertRaises(ValueError):
            fit(np.eye(3), np.zeros(3), np.ones(3))

    def test_negative_or_nonfinite_inputs_rejected(self):
        for y in [np.array([1, -1, 1]), np.array([1, np.nan, 1])]:
            with self.assertRaises(ValueError):
                fit(np.eye(3), y, np.ones(3))
        with self.assertRaises(ValueError):
            fit(np.eye(3), np.ones(3), np.zeros(3))

    @unittest.skipUnless((ROOT / 'data/raw/ref' / FILENAME).exists(), 'optional source XLSX not downloaded')
    def test_published_table_scale_and_exclusions(self):
        fine, cells, error, excluded = read_centroids(ROOT / 'data/raw/ref' / FILENAME)
        self.assertLess(error, 1e-8)
        self.assertEqual(cells.sum(), 20646)
        self.assertEqual(len(excluded), 27)
        self.assertEqual(fine.shape, (24448, 12))
        self.assertFalse(fine.index.has_duplicates)
        self.assertTrue(all(isinstance(gene, str) for gene in fine.index))
        # Date-corrupted rows are excluded AFTER validating full-table sums.
        self.assertTrue((fine.sum() < 1e6).all())


if __name__ == '__main__':
    unittest.main()
