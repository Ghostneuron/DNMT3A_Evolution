#!/usr/bin/env python3

from __future__ import annotations

import csv
import unittest
from pathlib import Path

from scripts.canonical_promoter_evolution import (
    CURATED_FULL_LENGTH_TSS_CLUSTER,
    PRIMARY_FULL_LENGTH,
    select_primary_targets,
    select_tss_cluster,
)


ROOT = Path(__file__).resolve().parents[1]


class CanonicalPromoterEvolutionTests(unittest.TestCase):
    def test_representatives_are_full_length_and_unique(self) -> None:
        with (
            ROOT / "results/06_isoform_evolution/target_transcript_tss.tsv"
        ).open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        selected = select_primary_targets(rows)
        self.assertEqual(len(selected), len(PRIMARY_FULL_LENGTH))
        self.assertTrue(all(
            row["contains_selected_full_length_CDS"] == "true"
            for row in selected
        ))
        self.assertEqual(
            {row["transcript_accession"] for row in selected},
            set(PRIMARY_FULL_LENGTH.values()),
        )
        cluster = select_tss_cluster(rows, selected)
        for species, accessions in CURATED_FULL_LENGTH_TSS_CLUSTER.items():
            self.assertEqual(
                {
                    row["transcript_accession"]
                    for row in cluster if row["species"] == species
                },
                accessions,
            )


if __name__ == "__main__":
    unittest.main()
