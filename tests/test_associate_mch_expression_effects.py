import unittest

from scripts.associate_mch_expression_effects import analyze


class MchExpressionAssociationTests(unittest.TestCase):
    def test_detects_inverse_effect_association(self):
        rows = []
        for chromosome in range(1, 5):
            for index in range(1, 11):
                mca = -index / 1000
                rows.append(
                    {
                        "chrom": f"chr{chromosome}",
                        "absolute_difference": str(mca),
                        "expression_difference": str(-2 * mca + chromosome / 10000),
                        "wt_mean_log2_CPM": str(index),
                        "wt_mean_corrected_mCA": str(index / 100),
                        "length_bp": str(1000 * index + chromosome),
                        "minimum_sites_across_samples": str(100 + index),
                    }
                )
        result = analyze(rows, bootstraps=20, seed=2)
        self.assertLess(result["unadjusted_spearman_rho"], -0.9)
        self.assertEqual(result["chromosome_blocks"], 4)


if __name__ == "__main__":
    unittest.main()
