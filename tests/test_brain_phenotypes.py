#!/usr/bin/env python3

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.import_brainzoo_mch import parse_brainzoo, phenotype_rows
from scripts.extract_evodevo_dnmt3a import stage_fields
from scripts.expression_trajectory_tests import (
    benjamini_hochberg,
    spearman,
)
from scripts.validate_brain_phenotypes import validate


class BrainPhenotypeTest(unittest.TestCase):
    def test_brainzoo_r_row_names_do_not_shift_columns(self) -> None:
        content = (
            "file\ttotal_CHcov\ttotal_mCH\tmCA\tmCT\tmCC\tspecies\n"
            "1\tsample.gz\t100\t4\t0.02\t0.01\t0.01\tMouse\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.tsv"
            path.write_text(content)
            rows = parse_brainzoo(path)
        self.assertEqual(rows[0]["file"], "sample.gz")
        self.assertEqual(rows[0]["species"], "Mouse")

    def test_mch_is_derived_from_coverage_counts(self) -> None:
        source = [{
            "file": "sample.gz", "total_CHcov": "100", "total_mCH": "4",
            "mCA": "0.02", "mCT": "0.01", "mCC": "0.01",
            "species": "Mouse",
        }]
        rows, audit = phenotype_rows(source, {"Mus musculus"})
        mch = next(row for row in rows if row["phenotype"].endswith("mCH_fraction"))
        self.assertEqual(float(mch["value"]), 0.04)
        self.assertEqual(audit[0]["panel_status"], "matched_panel")

    def test_validator_rejects_fraction_over_one(self) -> None:
        taxa = [{
            "species": "Mus musculus", "safe_id": "Mus_musculus",
            "ncbi_scientific_name": "Mus musculus", "is_mammal": "true",
        }]
        row = {
            "species": "Mus musculus", "phenotype": "global_brain_mCA_fraction",
            "value": "1.2", "unit": "fraction", "brain_region": "brain",
            "cell_type": "bulk tissue", "developmental_stage": "adult",
            "assay": "WGBS", "source_accession": "source",
            "citation_doi": "doi", "sample_id": "sample",
        }
        issues, _, summary = validate([row], taxa)
        self.assertTrue(any(i["code"] == "fraction_out_of_range" for i in issues))
        self.assertEqual(summary["errors"], 1)

    def test_stage_clocks_remain_distinct(self) -> None:
        self.assertEqual(stage_fields("e93"), ("embryo", "93", "day"))
        self.assertEqual(stage_fields("P23"), ("postnatal", "23", "day"))
        self.assertEqual(stage_fields("12wpb"), ("postnatal", "12", "week"))
        self.assertEqual(
            stage_fields("13.5", {"13.5": ("embryo", "day")}),
            ("embryo", "13.5", "day"),
        )
        self.assertEqual(
            stage_fields(
                "104", {"90": ("postnatal", "day")},
                postnatal_offset_days=14,
            ),
            ("postnatal", "90", "day"),
        )

    def test_spearman_and_bh(self) -> None:
        self.assertAlmostEqual(spearman([0, 1, 2], [3, 2, 1]), -1.0)
        adjusted = benjamini_hochberg([0.01, 0.04, 0.03])
        self.assertEqual(adjusted, [0.03, 0.04, 0.04])


if __name__ == "__main__":
    unittest.main()
