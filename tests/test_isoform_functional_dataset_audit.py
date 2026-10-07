import tempfile
import unittest
from pathlib import Path

from scripts.annotate_isoform_methylation_windows import (
    annotate_window,
    interval_distance,
)
from scripts.audit_isoform_functional_datasets import parse_soft_samples


class IsoformFunctionalDatasetAuditTests(unittest.TestCase):
    def test_soft_parser(self) -> None:
        text = """^SAMPLE = GSM1
!Sample_title = example
!Sample_source_name_ch1 = E15.5
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = genotype: Dnmt3a2-/-
!Sample_supplementary_file = https://example.org/a.gz
"""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "samples.soft"
            path.write_text(text)
            samples = parse_soft_samples(path)
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["sample_accession"], "GSM1")
        self.assertEqual(samples[0]["characteristics"]["tissue"], "brain")

    def test_interval_distance(self) -> None:
        self.assertEqual(interval_distance(10, 20, 15, 25), 0)
        self.assertEqual(interval_distance(10, 20, 20, 30), 0)
        self.assertEqual(interval_distance(10, 20, 25, 30), 5)

    def test_gene_annotation(self) -> None:
        genes = {"chr1": [(100, 200, "A"), (300, 400, "B")]}
        overlap, nearest, distance = annotate_window("chr1", 150, 250, genes)
        self.assertEqual(overlap, "A")
        self.assertEqual(nearest, "A")
        self.assertEqual(distance, 0)


if __name__ == "__main__":
    unittest.main()
