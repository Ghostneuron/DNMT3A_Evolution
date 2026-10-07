import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/04_selection/mammals"


class FinalCandidateEvidenceTests(unittest.TestCase):
    def setUp(self):
        if not (OUT / "final_candidate_evidence.json").exists():
            self.skipTest("Final candidate evidence has not been assembled")

    def test_candidate_ranking(self):
        summary = json.loads((OUT / "final_candidate_evidence.json").read_text())
        self.assertEqual(summary["candidate_order"][:3], ["G34", "T12", "S97"])
        self.assertEqual(
            summary["candidates_for_followup"], ["G34", "T12", "S97"]
        )
        self.assertEqual(
            set(summary["deprioritized_or_stopped"]), {"P5", "E18", "A114"}
        )

    def test_no_brain_claim_and_masked_statistics(self):
        with (OUT / "final_candidate_evidence.tsv").open() as handle:
            rows = {
                row["human_site"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertTrue(all(
            row["brain_phenotype_inference"] == "none" for row in rows.values()
        ))
        for site in ("G34", "T12"):
            self.assertGreaterEqual(
                float(rows[site]["masked_fubar_posterior"]), 0.9
            )
            self.assertLess(float(rows[site]["masked_meme_p"]), 0.05)
            self.assertGreater(float(rows[site]["masked_meme_adjusted_p"]), 0.05)
        self.assertEqual(
            rows["S97"]["annotation_alignment_effect"],
            "original_omission_caused_by_475_vertebrate_filter;"
            "100_percent_mammal_occupancy_and_four_alignment_agreement",
        )
        self.assertGreaterEqual(
            float(rows["S97"]["mammal_scope_macse_fubar_posterior"]), 0.9
        )
        self.assertAlmostEqual(
            float(rows["S97"]["mammal_scope_full_meme_p"]),
            0.02244551936189731,
        )
        self.assertGreater(
            float(rows["S97"]["mammal_scope_full_meme_bh_q_909"]), 0.05
        )

    def test_corrected_full_scan_does_not_override_reliability_gate(self):
        summary = json.loads((OUT / "final_candidate_evidence.json").read_text())
        self.assertEqual(
            summary["mammal_scope_full_meme_bh_discoveries"],
            ["P5", "E18", "S6"],
        )
        with (OUT / "final_candidate_evidence.tsv").open() as handle:
            rows = {
                row["human_site"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        for site in ("P5", "E18", "S6"):
            self.assertLess(
                float(rows[site]["mammal_scope_full_meme_bh_q_909"]), 0.05
            )
        self.assertEqual(rows["S6"]["recommended_action"], "exploratory_only")
        self.assertEqual(rows["P5"]["recommended_action"], "deprioritize")
        self.assertEqual(rows["E18"]["recommended_action"], "deprioritize")


if __name__ == "__main__":
    unittest.main()
