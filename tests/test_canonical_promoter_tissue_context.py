#!/usr/bin/env python3

from __future__ import annotations

import unittest

from scripts.fantom5_tissue_context import tissue_flags


class CanonicalPromoterTissueContextTests(unittest.TestCase):
    def test_brain_keyword_boundary_is_reused(self) -> None:
        self.assertEqual(tissue_flags("adult brain cortex"), (True, False))
        self.assertEqual(tissue_flags("fetal brain"), (True, True))
        self.assertEqual(tissue_flags("adrenal cortex"), (False, False))


if __name__ == "__main__":
    unittest.main()
