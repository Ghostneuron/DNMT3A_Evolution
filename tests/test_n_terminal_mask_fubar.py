import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/02_alignment/reliability_mask"


class NTerminalMaskFubarTests(unittest.TestCase):
    def setUp(self):
        if not (OUT / "fubar_summary.json").exists():
            self.skipTest("Masked FUBAR summary has not been generated")

    def test_retained_and_emergent_sites(self):
        summary = json.loads((OUT / "fubar_summary.json").read_text())
        self.assertEqual(set(summary["retained_primary_candidates"]), {"T12", "G34"})
        self.assertEqual(summary["mask_sensitivity_emergent"], ["A16"])

    def test_discovery_sites_do_not_reach_fubar_threshold(self):
        with (OUT / "candidate_fubar_comparison.tsv").open() as handle:
            rows = {
                row["human_site"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        for site in ("P5", "S6", "E18"):
            self.assertLess(float(rows[site]["masked_fubar_posterior"]), 0.9)
        self.assertGreaterEqual(float(rows["T12"]["masked_fubar_posterior"]), 0.9)
        self.assertGreaterEqual(float(rows["G34"]["masked_fubar_posterior"]), 0.9)


if __name__ == "__main__":
    unittest.main()
