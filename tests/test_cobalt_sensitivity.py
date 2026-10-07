import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/02_alignment/sensitivity/cobalt"


class CobaltSensitivityTests(unittest.TestCase):
    def test_prepared_input(self):
        path = OUT / "DNMT3A_mammals_unique_proteins.faa"
        if not path.exists():
            self.skipTest("COBALT input has not been prepared")
        self.assertEqual(path.read_text().count(">"), 250)

    def test_completed_alignment(self):
        summary_path = OUT / "COBALT_alignment_qc_summary.json"
        if not summary_path.exists():
            self.skipTest("COBALT alignment has not completed")
        summary = json.loads(summary_path.read_text())
        self.assertEqual(summary["sequences"], 250)
        self.assertEqual(summary["human_residues"], 912)
        self.assertEqual(summary["constraint_mode"], "norps_local_similarity_only")
        with (OUT / "COBALT_human_coordinate_crosswalk.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), summary["retained_human_mapped_columns"])
        self.assertTrue(all(int(row["human_DNMT3A1_aa"]) <= 912 for row in rows))

    def test_completed_selection_summary(self):
        path = ROOT / "results/04_selection/mammals/cobalt/cobalt_sensitivity_summary.json"
        if not path.exists():
            self.skipTest("COBALT selection sensitivity has not completed")
        summary = json.loads(path.read_text())
        self.assertEqual(
            summary["cobalt_fubar_candidates_ge_0.90"],
            ["T12", "G34", "S97"],
        )
        candidates = {
            row["human_site"]: row for row in summary["exact_meme_candidates"]
        }
        self.assertEqual(
            set(candidates),
            {"P5", "S6", "T12", "A16", "E18", "G34", "S97", "A114"},
        )
        self.assertLess(candidates["G34"]["holm_p_8"], 0.05)
        self.assertLess(candidates["T12"]["holm_p_8"], 0.05)


if __name__ == "__main__":
    unittest.main()
