import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from summarize_dunnart_neocortex_stages import (
    exact_permutation_p_greater,
    stage_from_sample,
    summarize_stages,
)


class SummarizeDunnartNeocortexStagesTests(unittest.TestCase):
    def test_stage_from_sample(self):
        self.assertEqual(stage_from_sample("P12 neocortex replicate 3"), "P12")

    def test_stage_pooling(self):
        records = [{
            "run_accession": "RUN1",
            "sample": "P12 neocortex replicate 1",
            "paired_read_count": 2_000_000,
            "junction_supporting_pair_counts": {
                "full_length_first_junction": 2,
                "internal_067_first_junction": 4,
                "internal_068_first_junction": 0,
                "common_downstream_junction": 8,
            },
        }]
        rows, summary = summarize_stages(records)
        self.assertEqual(rows[0]["internal_067_first_junction_pairs_per_million"], 2.0)
        self.assertEqual(
            summary["stages"]["P12"][
                "internal_067_fraction_of_diagnostic_first_junctions"
            ],
            0.6667,
        )
        self.assertEqual(
            summary["stages"]["P12"][
                "internal_067_to_common_downstream_pair_ratio"
            ],
            0.5,
        )
        self.assertEqual(
            summary["stages"]["P12"][
                "internal_067_to_full_length_first_pair_ratio"
            ],
            2.0,
        )

    def test_exact_permutation_complete_separation_three_per_group(self):
        self.assertEqual(
            exact_permutation_p_greater([4.0, 5.0, 6.0], [1.0, 2.0, 3.0]),
            0.05,
        )

    def test_exact_permutation_requires_balanced_groups(self):
        self.assertIsNone(exact_permutation_p_greater([2.0, 3.0], [1.0]))


if __name__ == "__main__":
    unittest.main()
