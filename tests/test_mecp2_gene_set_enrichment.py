import unittest

from scripts.analyze_mecp2_gene_set_enrichment import bh_adjust


class Mecp2GeneSetEnrichmentTests(unittest.TestCase):
    def test_bh_adjustment_is_monotone_in_rank(self):
        adjusted = bh_adjust([0.01, 0.04, 0.03, 0.5])
        self.assertAlmostEqual(adjusted[0], 0.04)
        self.assertAlmostEqual(adjusted[2], 0.05333333333333334)
        self.assertGreaterEqual(adjusted[3], adjusted[1])


if __name__ == "__main__":
    unittest.main()
