import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/07_functional_prioritization"


class FunctionalVariantPanelTests(unittest.TestCase):
    def test_tier_one_is_recurrent_and_lineage_supported(self):
        with (OUT / "DNMT3A1_variant_panel.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        tier_one = [row for row in rows if row["priority_tier"] == "1"]
        self.assertEqual(
            {row["experimental_substitution"] for row in tier_one},
            {"G34A", "G34S", "T12N", "T12A", "T12S", "S97G"},
        )
        self.assertTrue(all(
            int(row["change_event_count"]) >= 2
            and int(row["internal_branch_events"]) >= 1
            and row["brain_claim_allowed"] == "false"
            for row in tier_one
        ))

    def test_required_controls(self):
        with (OUT / "control_panel.tsv").open() as handle:
            controls = {
                row["control"] for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(
            controls,
            {"human_DNMT3A1_WT", "human_DNMT3A2", "matched_empty_vector"},
        )
        summary = json.loads((OUT / "summary.json").read_text())
        self.assertIn("neuronal", summary["interpretive_boundary"])

    def test_human_motor_cortex_context_is_calibration_only(self):
        summary = json.loads((OUT / "summary.json").read_text())
        context = summary["human_motor_cortex_assay_context"]
        self.assertEqual(context["reference_donors"], 2)
        self.assertGreater(
            context["neuronal_to_non_neuronal_median_ratio"], 3.0
        )
        with (OUT / "human_motor_cortex_assay_context.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(
            "candidate-variant effect" in row["interpretive_limit"]
            or "not independent biological replicates"
            in row["interpretive_limit"]
            for row in rows
        ))


if __name__ == "__main__":
    unittest.main()
