import unittest

import numpy as np

from scripts.model_mch_gene_body_covariates import analyze_contrast, ols, zscore


class MchGeneBodyCovariateModelTests(unittest.TestCase):
    def test_ols_recovers_linear_coefficients(self):
        x = np.arange(10, dtype=float)
        y = 2 + 3 * x
        beta, r_squared = ols(y, [x])
        self.assertAlmostEqual(beta[0], 2)
        self.assertAlmostEqual(beta[1], 3)
        self.assertAlmostEqual(r_squared, 1)

    def test_adjusted_length_effect_is_recovered(self):
        rows = []
        for chromosome in range(1, 5):
            for index in range(10):
                length = 1000 * (index + 1)
                wt = 0.01 + index * 0.001
                effect = -2 * wt - 0.002 * np.log10(length)
                rows.append(
                    {
                        "chrom": f"chr{chromosome}",
                        "length_bp": str(length),
                        "minimum_sites_across_samples": str(100 + index),
                        "wt_mean_corrected_mCA": str(wt),
                        "absolute_difference": str(effect),
                        "matched_WT_log2_CPM": str(index + 1),
                    }
                )
        result = analyze_contrast(rows, bootstraps=20, seed=1)
        self.assertLess(result["adjusted_standardized_length_beta"], 0)
        self.assertEqual(result["chromosome_blocks"], 4)

    def test_zscore_rejects_constant(self):
        with self.assertRaises(ValueError):
            zscore(np.ones(3))


if __name__ == "__main__":
    unittest.main()
