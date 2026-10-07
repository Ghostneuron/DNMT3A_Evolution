import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/02_alignment/reliability_mask"


class NTerminalMaskMemeTests(unittest.TestCase):
    def setUp(self):
        if not (OUT / "meme_summary.json").exists():
            self.skipTest("Masked MEME summary has not been generated")

    def test_nominal_but_not_corrected_results(self):
        summary = json.loads((OUT / "meme_summary.json").read_text())
        self.assertEqual(set(summary["nominal_p_le_0.05"]), {"T12", "G34"})
        self.assertEqual(summary["holm_significant_3_tests"], [])

    def test_a16_is_unsupported(self):
        with (OUT / "candidate_meme_comparison.tsv").open() as handle:
            rows = {
                row["human_site"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertGreater(float(rows["A16"]["meme_p"]), 0.1)
        self.assertEqual(rows["A16"]["classification"], "unsupported")


if __name__ == "__main__":
    unittest.main()
