import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/04_selection/mammals/s97_homology_audit"


class S97HomologyAuditTests(unittest.TestCase):
    def setUp(self):
        if not (OUT / "summary.json").exists():
            self.skipTest("S97 homology audit has not completed")

    def test_four_alignment_homology(self):
        summary = json.loads((OUT / "summary.json").read_text())
        self.assertEqual(summary["mammal_taxa"], 250)
        self.assertTrue(summary["all_pairwise_amino_acid_agreement"])
        self.assertEqual(
            set(summary["mammal_occupancy_all_alignments"].values()), {1.0}
        )

    def test_selection_is_consistent_but_nominal(self):
        with (OUT / "s97_selection_sensitivity.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["alignment"] for row in rows},
                         {"MACSE", "MAFFT", "PRANK", "COBALT"})
        self.assertTrue(all(
            float(row["fubar_positive_posterior"]) >= 0.9 for row in rows
        ))
        self.assertTrue(all(
            0.01 < float(row["meme_p"]) < 0.05 for row in rows
        ))

    def test_lineage_reconstruction(self):
        summary = json.loads((OUT / "lineages/summary.json").read_text())
        self.assertEqual(summary["root_residue"], "S")
        self.assertEqual(summary["terminal_reconstruction_mismatches"], 0)
        self.assertEqual(
            len(summary["prioritized_change_coupled_branches"]), 8
        )
        self.assertIn(
            "not branch-wide significance tests", summary["interpretation"]
        )


if __name__ == "__main__":
    unittest.main()
