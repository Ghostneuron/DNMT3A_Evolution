#!/usr/bin/env python3

from __future__ import annotations

import unittest

from scripts.summarize_esm_sensitivity import average_ranks, pearson, spearman


class EsmSensitivityTests(unittest.TestCase):
    def test_average_ranks_handles_ties(self) -> None:
        self.assertEqual(average_ranks([4, 1, 1, 8]), [3.0, 1.5, 1.5, 4.0])

    def test_correlations(self) -> None:
        self.assertAlmostEqual(pearson([1, 2, 3], [2, 4, 6]), 1.0)
        self.assertAlmostEqual(spearman([1, 3, 2], [10, 30, 20]), 1.0)


if __name__ == "__main__":
    unittest.main()
