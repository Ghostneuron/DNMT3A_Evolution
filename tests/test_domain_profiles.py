#!/usr/bin/env python3

from __future__ import annotations

import math
import unittest
from collections import Counter

from scripts.domain_profiles import hypergeometric_upper_tail, normalized_entropy


class DomainProfileMathTest(unittest.TestCase):
    def test_hypergeometric_upper_tail(self) -> None:
        observed = hypergeometric_upper_tail(10, 4, 2, 2)
        self.assertAlmostEqual(observed, math.comb(4, 2) / math.comb(10, 2))

    def test_normalized_entropy(self) -> None:
        self.assertEqual(normalized_entropy(Counter({"A": 10})), 0.0)
        expected = math.log(2) / math.log(20)
        self.assertAlmostEqual(normalized_entropy(Counter({"A": 5, "T": 5})), expected)


if __name__ == "__main__":
    unittest.main()
