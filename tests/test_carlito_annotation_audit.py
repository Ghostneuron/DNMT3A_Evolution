import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT
    / "results/04_selection/mammals/candidate_annotation_validation/Carlito_syrichta"
)


class CarlitoAnnotationAuditTests(unittest.TestCase):
    def setUp(self):
        if not (OUT / "summary.json").exists():
            self.skipTest("Carlito annotation audit has not been generated")

    def test_discordant_block(self):
        with (OUT / "discordant_site_exon_audit.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 22)
        self.assertEqual(
            (min(int(row["human_DNMT3A1_aa"]) for row in rows),
             max(int(row["human_DNMT3A1_aa"]) for row in rows)),
            (34, 59),
        )
        self.assertTrue(all(row["exon_number"] == "1" for row in rows))
        self.assertTrue(all(row["prank_codon"] == "---" for row in rows))
        self.assertTrue(all(
            row["decision"] == "exclude_from_homologous_site_inference"
            for row in rows
        ))

    def test_anchor_offset_and_model_evidence(self):
        summary = json.loads((OUT / "summary.json").read_text())
        self.assertEqual(summary["discordant_positions"], 22)
        self.assertEqual(summary["carlito_model_position_range"], [1, 22])
        self.assertEqual(summary["annotation"]["cds_start_in_rna_1based"], 97)
        self.assertEqual(summary["annotation"]["rna_genome_identity_pct"], 100.0)
        self.assertFalse(summary["annotation"]["independent_transcript_evidence"])


if __name__ == "__main__":
    unittest.main()
