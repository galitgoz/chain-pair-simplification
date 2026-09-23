"""Integrity checks for the prepared-input reproduction package."""
import csv
import hashlib
import json
import math
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


class PreparedInputTests(unittest.TestCase):
    def test_input_integrity(self):
        inputs = json.loads((ROOT/'data/inputs.json').read_text())
        self.assertEqual(len(inputs), 45)
        self.assertEqual(len({r['path'] for r in inputs}), 45)
        for row in inputs:
            with self.subTest(pair=row['pair']):
                self.assertEqual(hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest(), row['sha256'])

    def test_runnable_configs_match_reference(self):
        configs = json.loads((ROOT/'experiments/configurations.json').read_text())
        originals = {r['configuration_id']: r['records']['configuration.json']['content']
                     for r in json.loads((ROOT/'results/evidence.json').read_text())['attempts']}
        with (ROOT/'results/primary.csv').open(encoding='utf-8-sig', newline='') as f:
            rows = {r['configuration_id']: r for r in csv.DictReader(f) if r['observed'] == 'True'}
        self.assertEqual(len(configs), 180)
        self.assertEqual({c['configuration_id'] for c in configs}, set(rows))
        for c in configs:
            row = rows[c['configuration_id']]
            with self.subTest(config=c['configuration_id']):
                self.assertEqual(c['input_sha256'], row['input_sha256'])
                self.assertTrue((ROOT/c['input_file']).is_file())
                for key in ['delta1', 'delta2', 'delta3', 'alpha']:
                    self.assertEqual(c[key], originals[c['configuration_id']][key])
                    self.assertTrue(math.isclose(c[key], float(row[key]), rel_tol=1e-14, abs_tol=1e-12))
                self.assertFalse(c['certificate'])

    def test_consolidated_evidence_coverage(self):
        records = json.loads((ROOT/'results/evidence.json').read_text())['attempts']
        configs = json.loads((ROOT/'experiments/configurations.json').read_text())
        self.assertEqual(len(records), 180)
        self.assertEqual({r['configuration_id'] for r in records}, {c['configuration_id'] for c in configs})
        for record in records:
            self.assertIn('configuration.json', record['records'])
            self.assertIn('supervisor.json', record['records'])


if __name__ == '__main__':
    unittest.main()
