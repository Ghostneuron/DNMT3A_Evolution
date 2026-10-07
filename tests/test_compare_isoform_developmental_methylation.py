import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.compare_isoform_developmental_methylation import (
    analyze,
    build_contrasts,
    parse_beta,
)


class IsoformDevelopmentalMethylationTests(unittest.TestCase):
    def test_beta_parser(self) -> None:
        self.assertEqual(parse_beta("0.5"), 0.5)
        self.assertIsNone(parse_beta("NA"))
        with self.assertRaises(ValueError):
            parse_beta("1.2")

    def test_contrast_construction(self) -> None:
        headers = ["wt", "ko"]
        metadata = {
            "wt": {"stage": "E15.5", "tissue": "brain", "genotype": "WT"},
            "ko": {
                "stage": "E15.5",
                "tissue": "brain",
                "genotype": "Dnmt3a1-/-",
            },
        }
        contrasts = build_contrasts(headers, metadata)
        target = next(
            row for row in contrasts
            if row["contrast"] == "E15.5_brain_Dnmt3a1-/-_vs_WT"
        )
        self.assertEqual(target["WT_indexes"], [0])
        self.assertEqual(target["KO_indexes"], [1])

    def test_small_matrix(self) -> None:
        soft = """^SAMPLE = GSM1
!Sample_title = wt
!Sample_source_name_ch1 = E15.5
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = genotype: WT
^SAMPLE = GSM2
!Sample_title = ko
!Sample_source_name_ch1 = E15.5
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = genotype: Dnmt3a1-/-
"""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "samples.soft"
            metadata.write_text(soft)
            matrix = root / "beta.tsv.gz"
            with gzip.open(matrix, "wt") as handle:
                handle.write('"wt"\t"ko"\n')
                handle.write('"cg1"\t0.8\t0.4\n')
                handle.write('"cg2"\t0.6\t0.5\n')
            result = analyze(matrix=matrix, metadata_path=metadata)
        target = next(
            row for row in result["contrast_rows"]
            if row["contrast"] == "E15.5_brain_Dnmt3a1-/-_vs_WT"
        )
        self.assertAlmostEqual(target["mean_KO_minus_WT"], -0.25)


if __name__ == "__main__":
    unittest.main()
