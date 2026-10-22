from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from data_assets import safe_target


class AssetPathTests(unittest.TestCase):
    def test_project_relative_destination(self):
        result = safe_target({'path': 'data/raw/fixture.gz', 'url': 'https://example.org/fixture.gz'})
        self.assertEqual(result, ROOT / 'data/raw/fixture.gz')

    def test_absolute_and_parent_paths_rejected(self):
        for path in ['/tmp/fixture.gz', '../fixture.gz', 'data/../../fixture.gz', '.']:
            with self.assertRaises(ValueError):
                safe_target({'path': path, 'url': 'https://example.org/fixture.gz'})

    def test_non_https_source_rejected(self):
        with self.assertRaises(ValueError):
            safe_target({'path': 'data/raw/fixture.gz', 'url': 'file:///tmp/fixture.gz'})


if __name__ == '__main__':
    unittest.main()
