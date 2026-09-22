"""Layout checks requiring only the standard library."""
import importlib.util
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from repository_paths import resolve_recorded_path


class RepositoryLayoutTests(unittest.TestCase):
    def test_all_recorded_relocations_resolve(self):
        moves = json.loads((ROOT / 'docs/provenance/relocated_files.json').read_text())
        for old, new in moves.items():
            with self.subTest(path=old):
                self.assertFalse((ROOT / old).exists())
                self.assertTrue((ROOT / new).is_file())
                self.assertEqual(resolve_recorded_path(ROOT / old), ROOT / new)
                self.assertEqual(resolve_recorded_path(old), ROOT / new)

    def test_unrelated_paths_are_not_redirected(self):
        for path in [ROOT / 'README.md', ROOT / 'output/missing.py', ROOT.parent / 'curve_algorithms.py']:
            with self.subTest(path=path):
                self.assertEqual(resolve_recorded_path(path), path)

    def test_algorithm_import_locations(self):
        for name in ['curve_algorithms', 'cps_paper_algorithms', 'decimation_algorithms']:
            with self.subTest(module=name):
                self.assertEqual(Path(importlib.util.find_spec(name).origin), ROOT / 'src' / (name + '.py'))


if __name__ == '__main__':
    unittest.main()
