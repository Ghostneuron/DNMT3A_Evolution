import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.summarize_mch_contexts import (
    IndexedFasta,
    context_labels,
    oriented_trinucleotide,
    summarize,
)


class SummarizeMchContextsTests(unittest.TestCase):
    def test_orients_forward_and_reverse_cytosines(self):
        sequence = b"CAGACCGT"
        self.assertEqual(oriented_trinucleotide(sequence, 1), "CAG")
        self.assertEqual(oriented_trinucleotide(sequence, 7), "CGG")
        self.assertEqual(context_labels("CAG"), ("CA", "CHG", "CH", "all_C"))
        self.assertEqual(context_labels("CGG"), ("CG", "all_C"))

    def test_summarizes_context_counts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "test.fa"
            fasta.write_text(">chr1\nCAGCCATCG\n>phage_lambda\nCAAT\n")
            # name, length, offset, line bases, line width
            Path(f"{fasta}.fai").write_text(
                "chr1\t9\t6\t9\t10\nphage_lambda\t4\t30\t4\t5\n"
            )
            coverage = root / "calls.cov"
            coverage.write_text(
                "chr1\t1\t1\t50\t1\t1\n"
                "chr1\t5\t5\t25\t1\t3\n"
                "chr1\t8\t8\t100\t2\t0\n"
                "phage_lambda\t1\t1\t0\t0\t4\n"
            )
            genes = root / "genes.bed"
            genes.write_text("chr1\t0\t7\tGeneA\t0\t+\n")
            output = root / "out"
            report = summarize(
                coverage, fasta, output, bin_size=10, gene_bed=genes
            )
            self.assertEqual(report["input_rows"], 4)
            table = (output / "context_summary.tsv").read_text()
            self.assertIn("autosome\tCA\t2\t4\t2\t0.3333333333", table)
            self.assertIn("phage_lambda\tCA\t0\t4\t1\t0", table)
            with gzip.open(
                output / "gene_body_context_summary.tsv.gz", "rt"
            ) as handle:
                gene_table = handle.read()
            self.assertIn("GeneA\t+\tCA\t2\t4\t2\t0.3333333333", gene_table)

    def test_indexed_fasta_reuses_contig(self):
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / "x.fa"
            fasta.write_text(">x\nACGT\n")
            Path(f"{fasta}.fai").write_text("x\t4\t3\t4\t5\n")
            indexed = IndexedFasta(fasta)
            try:
                self.assertEqual(indexed.contig("x"), b"ACGT")
            finally:
                indexed.close()


if __name__ == "__main__":
    unittest.main()
