import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/04_selection/mammals/reliability_mask/lineages"


class MaskedCandidateLineageTests(unittest.TestCase):
    def setUp(self):
        if not (OUT / "summary.json").exists():
            self.skipTest("Masked candidate lineages have not been summarized")

    def test_summary_scope_and_reconstruction(self):
        summary = json.loads((OUT / "summary.json").read_text())
        sites = {row["human_site"]: row for row in summary["sites"]}
        self.assertEqual(set(sites), {"T12", "G34"})
        self.assertTrue(all(
            row["terminal_reconstruction_mismatches"] == 0
            for row in sites.values()
        ))
        self.assertGreaterEqual(sites["T12"]["masked_terminal_states_skipped"], 2)
        self.assertGreaterEqual(sites["G34"]["masked_terminal_states_skipped"], 1)

    def test_artifact_terminals_are_not_high_posterior(self):
        with (OUT / "prioritized_change_coupled_branches.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        observed = {(row["human_site"], row["branch"]) for row in rows}
        self.assertNotIn(("T12", "Puma_concolor"), observed)
        self.assertNotIn(("T12", "Vombatus_ursinus"), observed)
        self.assertNotIn(("G34", "Carlito_syrichta"), observed)


if __name__ == "__main__":
    unittest.main()
