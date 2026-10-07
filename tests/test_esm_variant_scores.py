#!/usr/bin/env python3

from __future__ import annotations

import unittest


class EsmVariantScoreSourceTests(unittest.TestCase):
    def test_script_documents_masked_marginal_boundary(self) -> None:
        from pathlib import Path

        script = (
            Path(__file__).resolve().parents[1]
            / "scripts/esm_variant_scores.py"
        ).read_text()
        self.assertIn("masked-marginal", script)
        self.assertIn("not a biochemical or adaptive-effect measurement", script)


if __name__ == "__main__":
    unittest.main()
