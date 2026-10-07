import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/04_selection/mammals/mammal_scope_filter_audit"


class MammalScopeFilterAuditTests(unittest.TestCase):
    def setUp(self):
        if not (OUT / "summary.json").exists():
            self.skipTest("Mammal-scope filter audit has not completed")

    def test_restored_sites_and_candidate(self):
        summary = json.loads((OUT / "summary.json").read_text())
        self.assertEqual(summary["restored_site_count"], 8)
        self.assertEqual(
            summary["restored_sites_fubar_ge_0_90_all_alignments"], ["S97"]
        )
        expected = {"T12", "G34", "S97"}
        for method in ("MACSE", "MAFFT", "COBALT"):
            self.assertEqual(
                set(summary["four_alignment_fubar_candidates"][method]),
                expected,
            )

    def test_other_restored_sites_are_not_candidates(self):
        with (OUT / "restored_human_sites.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(
            row["all_fubar_ge_0_90"] == (row["human_site"] == "S97").__str__().lower()
            for row in rows
        ))


if __name__ == "__main__":
    unittest.main()
