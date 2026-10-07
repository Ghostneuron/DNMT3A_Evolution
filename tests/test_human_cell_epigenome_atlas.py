#!/usr/bin/env python3

from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

from scripts.import_human_cell_epigenome_atlas import cell_class, quantile


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "results/05_brain_integration/human_cell_atlas_import_summary.json"


class HumanCellEpigenomeAtlasTest(unittest.TestCase):
    def test_cell_classification(self) -> None:
        self.assertEqual(
            cell_class("Neu Excitatory Intratelencephalic L5 1"),
            "neuron_excitatory",
        )
        self.assertEqual(
            cell_class("Glia Oligodendrocyte Progenitor"),
            "glia_oligodendrocyte_progenitor",
        )
        self.assertEqual(cell_class("Hema Myeloid Microglia"), "microglia")

    def test_linear_quantile(self) -> None:
        self.assertEqual(quantile([0.0, 1.0, 2.0, 3.0], 0.25), 0.75)
        self.assertEqual(quantile([0.0, 1.0, 2.0, 3.0], 0.75), 2.25)

    @unittest.skipUnless(
        (
            ROOT
            / "data/raw/human_cell_epigenome_atlas/5kCG100k3C_summary.csv.gz"
        ).exists(),
        "published atlas metadata not acquired",
    )
    def test_full_import_invariants(self) -> None:
        subprocess.run(
            [
                "/opt/anaconda3/bin/python",
                "scripts/import_human_cell_epigenome_atlas.py",
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        summary = json.loads(SUMMARY.read_text())
        self.assertEqual(summary["source_cells"], 86689)
        self.assertEqual(summary["tissues"], 16)
        self.assertEqual(summary["major_types"], 35)
        self.assertEqual(summary["annotated_subtypes"], 206)
        self.assertEqual(summary["unmapped_subtypes"], 0)
        self.assertEqual(summary["motor_cortex_cells"], 8225)
        self.assertGreater(
            summary["motor_cortex_neuronal_mCH_median"],
            summary["motor_cortex_non_neuronal_mCH_median"],
        )
        donor_ratios = summary[
            "motor_cortex_donor_neuronal_to_non_neuronal_mCH_median_ratio"
        ]
        self.assertEqual(set(donor_ratios), {"H1930001", "H1930002"})
        self.assertTrue(all(ratio > 3.0 for ratio in donor_ratios.values()))


if __name__ == "__main__":
    unittest.main()
