import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from dunnart_source_dnmt3a_expression import clean_header, cpm


class DunnartSourceDnmt3aExpressionTests(unittest.TestCase):
    def test_clean_header(self):
        self.assertEqual(clean_header('"FTD_P12_436.genes.results"'), "FTD_P12_436.genes.results")

    def test_cpm(self):
        self.assertEqual(cpm(25, 2_500_000), 10)


if __name__ == "__main__":
    unittest.main()
