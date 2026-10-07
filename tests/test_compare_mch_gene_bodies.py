import csv
import gzip
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.compare_mch_gene_bodies import build, corrected


SAMPLES = {
    "w1": {"genotype": "WT"},
    "w2": {"genotype": "WT"},
    "a1": {"genotype": "Dnmt3a1_KO"},
    "a2": {"genotype": "Dnmt3a1_KO"},
    "b1": {"genotype": "Dnmt3a2_KO"},
    "b2": {"genotype": "Dnmt3a2_KO"},
    "d1": {"genotype": "Dnmt3a1_delta_N"},
    "d2": {"genotype": "Dnmt3a1_delta_N"},
}


class GeneBodyComparisonTests(unittest.TestCase):
    def test_background_correction(self):
        self.assertAlmostEqual(corrected(0.02, 0.01), 0.0101010101)
        self.assertEqual(corrected(0.005, 0.01), 0.0)

    def test_non_autosomal_gene_is_excluded(self):
        # The primary endpoint and downstream chromosome bootstrap are chr1--chr19.
        from scripts.compare_mch_gene_bodies import AUTOSOMES

        self.assertIn("chr19", AUTOSOMES)
        self.assertNotIn("chrX", AUTOSOMES)

    def test_complete_case_gene_contrasts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fractions = {"w1": .02, "w2": .022, "a1": .01, "a2": .01,
                         "b1": .021, "b2": .023, "d1": .015, "d2": .016}
            for run, fraction in fractions.items():
                source = root / "results/bismark" / run / "context_summary"
                source.mkdir(parents=True)
                (source / "context_summary.tsv").write_text(
                    "scope\tcontext\tmethylated\tunmethylated\tsites\tfraction\n"
                    "phage_lambda\tCA\t10\t990\t250\t0.01\n"
                )
                with gzip.open(source / "gene_body_context_summary.tsv.gz", "wt", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(["chrom", "start", "end", "gene", "strand", "context", "methylated", "unmethylated", "sites", "fraction"])
                    writer.writerow(["chr1", 0, 200000, "GeneA", "+", "CA", 20, 980, 250, fraction])
            output = root / "out"
            with patch("scripts.compare_mch_gene_bodies.read_samples", return_value=SAMPLES):
                report = build(root, output, min_sites=100)
            self.assertTrue(report["complete"])
            self.assertEqual(report["complete_case_genes"], 1)
            text = (output / "gene_body_mCA_summary.tsv").read_text()
            self.assertIn("Dnmt3a1_KO-WT", text)
            with (output / "gene_body_mCA_summary.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["genes"], "1")
            self.assertEqual(rows[0]["genes_ge_100kb"], "1")


if __name__ == "__main__":
    unittest.main()
