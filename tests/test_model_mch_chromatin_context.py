import unittest

import numpy as np

from scripts.model_mch_chromatin_context import add_enrichment_features, analyze_feature


class MchChromatinContextModelTests(unittest.TestCase):
    def test_enrichment_is_log1p_chip_minus_input(self):
        row = {}
        for region in ("gene_body", "promoter"):
            row[f"input_{region}_mean"] = "1"
            row[f"dnmt3a_flag_{region}_mean"] = "3"
            row[f"h3k4me3_{region}_mean"] = "1"
            row[f"h3k27me3_{region}_mean"] = "0"
        values = add_enrichment_features(row)
        self.assertAlmostEqual(values["dnmt3a_flag_gene_body_enrichment"], np.log(2))
        self.assertAlmostEqual(values["h3k4me3_promoter_enrichment"], 0)

    def test_adjusted_focal_signal_is_recovered(self):
        rows = []
        for chromosome in range(1, 6):
            for index in range(12):
                focal = (-1 if index % 2 else 1) * (1 + chromosome * 0.03)
                wt = 0.004 + index * 0.0007
                length = 1000 * (index + 2)
                rows.append(
                    {
                        "chrom": f"chr{chromosome}",
                        "length_bp": length,
                        "minimum_sites_across_samples": 100 + index * 7,
                        "wt_mean_corrected_mCA": wt,
                        "matched_WT_log2_CPM": index + chromosome / 10,
                        "absolute_difference": -wt + 0.004 * focal,
                        "dnmt3a_flag_gene_body_enrichment": focal,
                    }
                )
        result = analyze_feature(
            rows, "dnmt3a_flag_gene_body_enrichment", bootstraps=30, seed=7
        )
        self.assertGreater(result["standardized_feature_beta"], 0)
        self.assertGreater(result["feature_partial_r_squared"], 0)
        self.assertEqual(result["chromosome_blocks"], 5)


if __name__ == "__main__":
    unittest.main()
