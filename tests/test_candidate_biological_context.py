import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLE = (
    ROOT
    / "results/04_selection/mammals/candidate_biological_context"
    / "retained_candidate_context.tsv"
)


class CandidateBiologicalContextTests(unittest.TestCase):
    def setUp(self):
        if not TABLE.exists():
            self.skipTest("Candidate biological context has not been built")

    def test_retained_candidates_are_dnmt3a1_disordered_tail_sites(self):
        with TABLE.open() as handle:
            rows = {
                row["human_site"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(set(rows), {"G34", "T12", "S97"})
        self.assertTrue(all(
            row["uniprot_disordered_region"] == "true"
            and row["dnmt3a1_specific_absent_from_isoform2"] == "true"
            and float(row["alphafold_plddt"]) < 50
            for row in rows.values()
        ))
        self.assertEqual(rows["S97"]["exact_uniprot_modified_residue"], "")
        self.assertEqual(
            rows["S97"]["nearest_uniprot_modified_residue"],
            "Phosphoserine@105",
        )


if __name__ == "__main__":
    unittest.main()
