import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from dunnart_stage_robustness import contrast, pooled_metrics


def row(stage, pairs, junction_internal, common, junction_full, salmon_input, retained,
        salmon_internal, internal_full_ratio):
    return {
        "stage": stage,
        "paired_read_count": pairs,
        "internal_067_first_junction_pairs": junction_internal,
        "common_downstream_junction_pairs": common,
        "full_length_first_junction_pairs": junction_full,
        "input_fragments": salmon_input,
        "retained_dnmt3a_fragments": retained,
        "internal_066_067_estimated_fragments": salmon_internal,
        "internal_066_067_to_full_length_count_ratio": internal_full_ratio,
    }


class DunnartStageRobustnessTests(unittest.TestCase):
    def test_pooled_metrics_use_pooled_denominators(self):
        rows = [
            row("P12", 100, 10, 20, 5, 100, 20, 10, 2),
            row("P12", 300, 15, 30, 5, 300, 30, 15, 3),
        ]
        metrics = pooled_metrics(rows)
        self.assertEqual(metrics["junction_internal_rate"], 25 / 400)
        self.assertEqual(metrics["junction_internal_common_ratio"], 25 / 50)
        self.assertEqual(metrics["junction_internal_full_ratio"], 25 / 10)
        self.assertEqual(metrics["salmon_internal_rate"], 25 / 400)
        self.assertEqual(metrics["salmon_internal_fraction"], 25 / 50)
        self.assertEqual(metrics["salmon_internal_full_ratio"], 2.5)

    def test_contrast_is_p12_over_p20(self):
        rows = [
            row("P12", 100, 20, 40, 10, 100, 20, 10, 2),
            row("P20", 100, 10, 40, 10, 100, 20, 5, 2),
        ]
        values = contrast(rows)
        self.assertEqual(values["junction_internal_rate_fold"], 2)
        self.assertEqual(values["junction_internal_common_ratio_fold"], 2)
        self.assertEqual(values["junction_internal_full_ratio_fold"], 2)
        self.assertEqual(values["salmon_internal_rate_fold"], 2)
        self.assertEqual(values["salmon_internal_fraction_fold"], 2)
        self.assertEqual(values["salmon_internal_full_ratio_fold"], 1)


if __name__ == "__main__":
    unittest.main()
