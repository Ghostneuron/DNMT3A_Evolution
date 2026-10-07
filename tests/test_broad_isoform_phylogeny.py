#!/usr/bin/env python3

from __future__ import annotations

import unittest

from scripts.broad_isoform_phylogeny import broad_present, strict_present


class BroadIsoformPhylogenyTests(unittest.TestCase):
    def test_presence_rules_do_not_treat_absence_as_loss(self) -> None:
        row = {
            "A2_like_products": "0",
            "broader_downstream_core_status": (
                "no_downstream_core_product_detected"
            ),
        }
        self.assertFalse(strict_present(row))
        self.assertFalse(broad_present(row))
        row["broader_downstream_core_status"] = "leaderless_core_present"
        self.assertTrue(broad_present(row))
        self.assertFalse(strict_present(row))


if __name__ == "__main__":
    unittest.main()
