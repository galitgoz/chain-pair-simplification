"""Regression checks for saved-result selection, not optimizer correctness."""
import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import summarize_primary as summary


class SavedResultChecks(unittest.TestCase):
    def test_recorded_paper_results(self):
        rows, medians, result = summary.summarize()
        self.assertEqual(len(rows), 172)
        self.assertEqual(result['quality'], {
            'equal_K_IC_CPS2F': 38, 'CPS2F_K_plus_one': 5,
            'IC_coupling_violations': 29, 'IC_violations_repaired_at_equal_K': 24,
            'CPS2F_smaller_K_than_CPS3F': 42, 'CPS2F_coupling_satisfied': 43})
        self.assertEqual({(r['domain'], r['method']): round(r['median_solver_seconds'], 3)
                          for r in medians}, {
            ('protein', 'CPS-2F'): 1.858, ('protein', 'CPS-3F'): 0.032,
            ('hurricane', 'CPS-2F'): 7.442, ('hurricane', 'CPS-3F'): 0.039})

    def altered_source(self, mutation, message):
        with summary.SOURCE.open(encoding='utf-8-sig', newline='') as stream:
            rows = list(csv.DictReader(stream))
        mutation(rows)
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'results.csv'
            summary.write_csv(source, rows)
            with patch.object(summary, 'SOURCE', source):
                with self.assertRaisesRegex(ValueError, message):
                    summary.summarize()

    def test_duplicate_primary_observation_is_rejected(self):
        self.altered_source(lambda rows: rows.append(next(r.copy() for r in rows if r['validated'] == 'True')),
                            'duplicate method')

    def test_mismatched_threshold_is_rejected(self):
        def mutate(rows):
            row = next(r for r in rows if r['validated'] == 'True' and r['domain'] == 'protein')
            row['delta1'] = str(float(row['delta1']) + 1)
        self.altered_source(mutate, 'do not share inputs and thresholds')


if __name__ == '__main__':
    unittest.main()
