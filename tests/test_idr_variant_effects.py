#!/usr/bin/env python3

from __future__ import annotations

import unittest

from scripts.idr_variant_effects import (
    bh_adjust,
    mutate,
    normalized_entropy,
    sequence_charge_decoration,
    sequence_metrics,
    window,
)


class IdrVariantEffectMathTests(unittest.TestCase):
    def test_bh_adjustment(self) -> None:
        adjusted = bh_adjust([0.01, 0.04, 0.03])
        self.assertAlmostEqual(adjusted[0], 0.03)
        self.assertAlmostEqual(adjusted[1], 0.04)
        self.assertAlmostEqual(adjusted[2], 0.04)

    def test_mutate_requires_matching_reference(self) -> None:
        self.assertEqual(mutate("ABCDE", 3, "C", "G"), "ABGDE")
        with self.assertRaises(RuntimeError):
            mutate("ABCDE", 3, "A", "G")

    def test_terminal_window_is_clipped(self) -> None:
        start, end, local = window("ABCDEFGHIJ", 1, 5)
        self.assertEqual((start, end, local), (1, 3, "ABC"))

    def test_entropy_and_charge_decoration(self) -> None:
        self.assertEqual(normalized_entropy("AAAA"), 0.0)
        self.assertEqual(sequence_charge_decoration("AAAA"), 0.0)

    def test_metrics_detect_phospho_acceptor_loss(self) -> None:
        reference = sequence_metrics("ASTP")
        alternative = sequence_metrics("ASAP")
        self.assertGreater(
            reference["phospho_acceptor_fraction"],
            alternative["phospho_acceptor_fraction"],
        )


if __name__ == "__main__":
    unittest.main()
