import unittest

import numpy as np

from scripts.model_p21_polycomb_context import analyze


class P21PolycombModelTests(unittest.TestCase):
    def test_analyze_recovers_adjusted_feature(self):
        chromosomes = np.repeat([f"chr{i}" for i in range(1, 6)], 12)
        index = np.arange(len(chromosomes), dtype=float)
        base = index / len(index)
        focal = np.where((index.astype(int) % 2) == 0, 1.0, -1.0)
        outcome = 2 * base - 0.5 * focal
        result = analyze(outcome, [base], [focal], ["mark"], chromosomes, 30, 5)
        self.assertLess(result["mark_standardized_beta"], 0)
        self.assertGreater(result["joint_partial_r_squared"], 0)
        self.assertEqual(result["chromosome_blocks"], 5)


if __name__ == "__main__":
    unittest.main()
