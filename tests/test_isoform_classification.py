#!/usr/bin/env python3

from __future__ import annotations

import unittest

from scripts.classify_dnmt3a_isoforms import classify_product


class IsoformClassificationTest(unittest.TestCase):
    def test_operational_a2_junction(self) -> None:
        full = "A" * 213 + "C" * 699
        a2 = "M" * 24 + "C" * 699
        product_class, _ = classify_product(full, a2)
        self.assertEqual(product_class, "DNMT3A2_like_downstream_start")

    def test_exact_suffix_without_unique_leader_is_not_a2(self) -> None:
        full = "A" * 223 + "C" * 689
        suffix = "C" * 689
        product_class, _ = classify_product(full, suffix)
        self.assertEqual(product_class, "leaderless_downstream_core_product")


if __name__ == "__main__":
    unittest.main()
