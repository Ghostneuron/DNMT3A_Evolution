#!/usr/bin/env python3

from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "results/08_synthesis/project_readiness.json"


class ProjectReadinessAuditTest(unittest.TestCase):
    def test_claim_boundaries_and_next_actions(self):
        subprocess.run(
            ["/opt/anaconda3/bin/python", "scripts/project_readiness_audit.py"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        summary = json.loads(SUMMARY.read_text())
        self.assertEqual(summary["computational_discovery_phase"], "complete")
        self.assertFalse(summary["brain_association_ready"])
        self.assertEqual(
            summary["coding_candidates_for_followup"], ["G34", "T12", "S97"]
        )
        self.assertEqual(len(summary["tier_1_functional_substitutions"]), 6)
        brain = next(
            row for row in summary["readiness_rows"]
            if row["inference_layer"] == "comparative_brain_association"
        )
        self.assertEqual(brain["status"], "blocked_on_new_species_data")
        self.assertIn("brain adaptation", brain["not_supported"])
        functional = next(
            row for row in summary["readiness_rows"]
            if row["inference_layer"] == "functional_prioritization"
        )
        self.assertEqual(
            functional["status"],
            "in_silico_screen_complete_no_exceptional_candidates",
        )
        self.assertGreater(
            summary["in_silico_functional_screen"][
                "esm_natural_null_spearman_rho"
            ],
            0.95,
        )
        promoter = next(
            row for row in summary["readiness_rows"]
            if row["inference_layer"] == "dnmt3a1_promoter_architecture"
        )
        self.assertEqual(
            promoter["status"], "human_mouse_cluster_orthology_supported"
        )
        self.assertEqual(
            summary["dnmt3a1_promoter_screen"][
                "canonical_tss_cluster_models"
            ],
            8,
        )
        lineage = next(
            row for row in summary["readiness_rows"]
            if row["inference_layer"]
            == "lineage_specific_promoter_isoform_evolution"
        )
        self.assertEqual(
            lineage["status"],
            "therian_ancestry_supported_human_acceleration_not_detected",
        )
        self.assertEqual(
            summary["dense_promoter_isoform_evolution"][
                "eutherian_orders_with_any_downstream_core_product"
            ],
            13,
        )
        self.assertIn("post hoc selection tests", summary["stop_rule"])


if __name__ == "__main__":
    unittest.main()
