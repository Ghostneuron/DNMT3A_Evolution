import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fantom5_tissue_context import decode_sample, tissue_flags


class Fantom5TissueContextTests(unittest.TestCase):
    def test_sample_decoding(self):
        accession, description = decode_sample(
            "counts.fetal%20brain%2c%20donor1.CNhs12345.library.hg38.nobarcode"
        )
        self.assertEqual(accession, "CNhs12345")
        self.assertEqual(description, "fetal brain, donor1")

    def test_brain_and_development_keywords(self):
        self.assertEqual(tissue_flags("fetal brain, donor1"), (True, True))
        self.assertEqual(tissue_flags("cortex, neonate N30"), (True, True))
        self.assertEqual(
            tissue_flags("Renal Cortical Epithelial Cells, donor1"),
            (False, False),
        )
        self.assertEqual(tissue_flags("adult liver, donor1"), (False, False))


if __name__ == "__main__":
    unittest.main()
