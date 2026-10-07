import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/04_selection/mammals/candidate_annotation_validation"


class CandidateAnnotationAuditTests(unittest.TestCase):
    def setUp(self):
        if not (OUT / "candidate_annotation_summary.json").exists():
            self.skipTest("Candidate annotation audit has not been run")

    def test_exon_mapping_and_decisions(self):
        with (OUT / "candidate_site_exon_audit.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 8)
        mapped = [row for row in rows if row["model_residue"] != "-"]
        self.assertTrue(all(row["crosses_exon_junction"] == "no" for row in mapped))
        puma_early = [
            row for row in rows
            if row["species"] == "Puma_concolor"
            and row["human_numbered_site"] in {"P5", "S6", "E18"}
        ]
        self.assertTrue(all(row["exon_number"] == "1" for row in puma_early))
        self.assertTrue(all(
            row["positional_homology_decision"]
            == "exclude_from_homologous_site_inference"
            for row in puma_early
        ))
        self.assertEqual(
            [(row["model_aa_position"], row["model_residue"]) for row in puma_early],
            [("2", "R"), ("3", "P"), ("15", "E")],
        )
        puma_g34 = next(
            row for row in rows
            if row["species"] == "Puma_concolor"
            and row["human_numbered_site"] == "G34"
        )
        self.assertEqual(puma_g34["model_residue"], "-")
        self.assertEqual(
            puma_g34["positional_homology_decision"], "alignment_gap_no_residue"
        )

    def test_annotation_model_summary(self):
        with (OUT / "annotation_model_audit.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            {row["species"]: row["exon_count"] for row in rows},
            {"Puma_concolor": "23", "Vombatus_ursinus": "22"},
        )
        self.assertTrue(all(row["rna_genome_identity_pct"] == "100.0" for row in rows))
        self.assertTrue(all("GC:1" in row["donor_sites"] for row in rows))

    def test_puma_relative_anchor(self):
        with (OUT / "n_terminal_relative_anchor_audit.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        puma = [row for row in rows if row["target_species"] == "Puma_concolor"]
        self.assertEqual(len(puma), 4)
        self.assertTrue(all(row["target_match_start_1based"] == "26" for row in puma))
        self.assertTrue(all(row["relative_match_start_1based"] == "57" for row in puma))
        self.assertTrue(all(int(row["exact_match_length_aa"]) >= 124 for row in puma))

    def test_summary_preserves_inference_boundary(self):
        summary = json.loads((OUT / "candidate_annotation_summary.json").read_text())
        self.assertIn("not orthology", summary["inference_boundary"])
        self.assertEqual(
            summary["site_decision_counts"]["exclude_from_homologous_site_inference"],
            3,
        )


if __name__ == "__main__":
    unittest.main()
