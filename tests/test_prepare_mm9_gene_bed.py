import tempfile
import unittest
from pathlib import Path

from scripts.prepare_mm9_gene_bed import prepare


class PrepareMm9GeneBedTests(unittest.TestCase):
    def test_selects_longest_protein_coding_transcript(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "refGene.txt"
            output = root / "genes.bed"
            source.write_text(
                "0\tNM_short\tchr1\t+\t100\t200\t100\t200\t1\t100,\t200,\t0\tGeneA\n"
                "0\tNM_long\tchr1\t+\t50\t300\t50\t300\t1\t50,\t300,\t0\tGeneA\n"
                "0\tNR_noncoding\tchr1\t+\t0\t500\t0\t0\t1\t0,\t500,\t0\tGeneB\n"
            )
            count = prepare(source, output)
            self.assertEqual(count, 1)
            self.assertEqual(output.read_text(), "chr1\t50\t300\tGeneA\t0\t+\n")


if __name__ == "__main__":
    unittest.main()
