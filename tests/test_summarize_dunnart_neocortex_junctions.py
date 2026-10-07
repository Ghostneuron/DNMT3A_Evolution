import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from summarize_dunnart_neocortex_junctions import summarize


class SummarizeDunnartNeocortexJunctionTests(unittest.TestCase):
    def test_pooled_counts_and_rates(self):
        records = [
            {
                "run_accession": "RUN1",
                "sample": "rep1",
                "paired_read_count": 1_000_000,
                "junction_supporting_pair_counts": {
                    "full_length_first_junction": 1,
                    "internal_067_first_junction": 3,
                    "internal_068_first_junction": 0,
                    "common_downstream_junction": 4,
                },
            },
            {
                "run_accession": "RUN2",
                "sample": "rep2",
                "paired_read_count": 1_000_000,
                "junction_supporting_pair_counts": {
                    "full_length_first_junction": 1,
                    "internal_067_first_junction": 1,
                    "internal_068_first_junction": 0,
                    "common_downstream_junction": 2,
                },
            },
        ]
        rows, pooled = summarize(records)
        self.assertEqual(len(rows), 2)
        self.assertEqual(
            pooled["pooled_junction_supporting_pair_counts"][
                "internal_067_first_junction"
            ],
            4,
        )
        self.assertEqual(
            pooled["pooled_pairs_per_million"][
                "internal_067_first_junction"
            ],
            2.0,
        )
        self.assertEqual(
            pooled["internal_067_fraction_of_diagnostic_first_junctions"],
            0.6667,
        )


if __name__ == "__main__":
    unittest.main()
