#!/usr/bin/env python3

from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = (
    ROOT
    / "results/05_brain_integration/human_mch_cross_dataset_calibration.json"
)


class HumanMchCalibrationTest(unittest.TestCase):
    def test_independent_bulk_value_is_bracketed_without_mixture_claim(self):
        subprocess.run(
            ["/opt/anaconda3/bin/python", "scripts/calibrate_human_mch_context.py"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        summary = json.loads(SUMMARY.read_text())
        self.assertTrue(summary["bulk_value_bracketed_by_atlas_cell_classes"])
        self.assertGreater(
            summary["brainzoo_bulk_mCH_fraction"],
            summary["atlas_non_neuronal_mCH_median"],
        )
        self.assertLess(
            summary["brainzoo_bulk_mCH_fraction"],
            summary["atlas_neuronal_mCH_median"],
        )
        self.assertIn("Do not estimate neuronal fraction", summary["forbidden_inference"])


if __name__ == "__main__":
    unittest.main()
