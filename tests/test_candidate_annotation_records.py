import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/04_selection/mammals/candidate_annotation_validation"


class CandidateAnnotationRecordTests(unittest.TestCase):
    def test_manifest(self):
        path = OUT / "record_manifest.tsv"
        if not path.exists():
            self.skipTest("Candidate annotation records have not been extracted")
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            {row["safe_id"] for row in rows},
            {"Carlito_syrichta", "Puma_concolor", "Vombatus_ursinus"},
        )
        self.assertTrue(all(row["transcript_type"] == "PROTEIN_CODING_MODEL" for row in rows))
        self.assertTrue(all(row["transcript_count"] == "1" for row in rows))
        self.assertTrue(all(row["scaffold_type"] == "Unplaced Scaffold" for row in rows))


if __name__ == "__main__":
    unittest.main()
