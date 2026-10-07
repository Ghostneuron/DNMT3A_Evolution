import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/05_brain_integration"


class CandidateBrainPhenotypeCoverageTests(unittest.TestCase):
    def test_no_association_is_overcalled(self):
        summary = json.loads(
            (OUT / "candidate_brain_phenotype_coverage_summary.json").read_text()
        )
        self.assertEqual(summary["candidate_sites"], ["T12", "G34", "S97"])
        self.assertEqual(summary["association_ready_combinations"], [])
        self.assertEqual(
            summary["s97_nonreference_taxa_with_any_brain_phenotype"],
            ["Ornithorhynchus_anatinus"],
        )

    def test_expression_overlay_includes_s97(self):
        summary = json.loads(
            (OUT / "candidate_expression_overlay_summary.json").read_text()
        )
        self.assertEqual(summary["S97_state_counts"], {"S": 6})
        self.assertFalse(summary["candidate_expression_association_testable"])
        with (OUT / "candidate_expression_overlay.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(all(
            row["mammal_scope_MACSE_S97"] == row["PRANK_S97"] == "S"
            for row in rows
        ))


if __name__ == "__main__":
    unittest.main()
