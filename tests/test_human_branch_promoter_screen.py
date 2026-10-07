#!/usr/bin/env python3

from __future__ import annotations

import unittest

from scripts.human_branch_promoter_screen import parse_maf_sequences


class HumanBranchPromoterScreenTests(unittest.TestCase):
    def test_semicolon_encoded_maf_parser(self) -> None:
        block = (
            "a score=1;"
            "s hg38.chr2 10 3 + 100 ACG;"
            "s panTro6.chr2A 12 3 + 100 ATG;"
            "i panTro6.chr2A C 0 C 0"
        )
        parsed = parse_maf_sequences(block)
        self.assertEqual(parsed["hg38"][1], "ACG")
        self.assertEqual(parsed["panTro6"][1], "ATG")


if __name__ == "__main__":
    unittest.main()
