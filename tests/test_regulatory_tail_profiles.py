#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import unittest
from pathlib import Path

from scripts.regulatory_tail_profiles import (
    hypergeometric_lower_tail,
    odds_ratio,
)


class RegulatoryTailProfileMathTest(unittest.TestCase):
    def test_hypergeometric_lower_tail(self) -> None:
        observed = hypergeometric_lower_tail(10, 4, 2, 0)
        self.assertAlmostEqual(observed, math.comb(6, 2) / math.comb(10, 2))

    def test_odds_ratio(self) -> None:
        self.assertAlmostEqual(odds_ratio(2, 8, 4, 6), 0.375)


class RegulatoryTailProfileOutputTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        root = Path(__file__).resolve().parents[1]
        output = (
            root
            / "results/04_selection/mammals/regulatory_tail/"
            "regulatory_subregion_summary.tsv"
        )
        with output.open() as handle:
            cls.rows = {
                row["region"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }

    def test_engagement_region_is_conserved_relative_to_upstream(self) -> None:
        engagement = self.rows["nucleosome_H2AK119ub_engagement"]
        upstream = self.rows["upstream_disordered_tail"]
        self.assertEqual(int(engagement["retained_sites"]), 56)
        self.assertEqual(int(engagement["polymorphic_sites"]), 6)
        self.assertEqual(int(upstream["retained_sites"]), 154)
        self.assertEqual(int(upstream["polymorphic_sites"]), 132)
        self.assertLess(
            float(engagement["mean_normalized_amino_acid_entropy"]),
            float(upstream["mean_normalized_amino_acid_entropy"]),
        )


if __name__ == "__main__":
    unittest.main()
