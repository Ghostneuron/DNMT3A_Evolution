#!/usr/bin/env python3
"""Summarize FANTOM5 tissue detection of canonical DNMT3A1 TSS clusters."""

from __future__ import annotations

import csv
import json
from pathlib import Path

try:
    from scripts.fantom5_tissue_context import (
        MATRICES,
        decode_sample,
        detected,
        read_target_counts,
        tissue_flags,
    )
except ModuleNotFoundError:
    from fantom5_tissue_context import (
        MATRICES,
        decode_sample,
        detected,
        read_target_counts,
        tissue_flags,
    )


ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT / "results/06_isoform_evolution/canonical_dnmt3a1_promoter"
)
NEARBY = OUT / "canonical_fantom5_nearby_cage_peaks.tsv"


def main() -> None:
    with NEARBY.open() as handle:
        nearby = list(csv.DictReader(handle, delimiter="\t"))
    sample_rows = []
    summaries = []
    for species, matrix in MATRICES.items():
        local_peaks = {
            row["peak_id"] for row in nearby if row["species"] == species
        }
        columns, counts = read_target_counts(matrix, local_peaks)
        for index, column in enumerate(columns):
            accession, description = decode_sample(column)
            brain, developmental = tissue_flags(description)
            cluster_count = sum(values[index] for values in counts.values())
            sample_rows.append({
                "species": species,
                "sample_accession": accession,
                "sample_description": description,
                "brain_CNS_keyword_match": str(brain).lower(),
                "developmental_keyword_match": str(developmental).lower(),
                "canonical_cluster_peak_count": len(local_peaks),
                "canonical_cluster_raw_count": cluster_count,
                "canonical_cluster_detected": str(cluster_count > 0).lower(),
                "quantification_limit": (
                    "Raw counts used for within-library detection; magnitudes "
                    "are not compared across unnormalized libraries"
                ),
            })
        species_rows = [
            row for row in sample_rows if row["species"] == species
        ]
        brain_rows = [
            row for row in species_rows
            if row["brain_CNS_keyword_match"] == "true"
        ]
        nonbrain_rows = [
            row for row in species_rows
            if row["brain_CNS_keyword_match"] == "false"
        ]
        developmental_brain = [
            row for row in brain_rows
            if row["developmental_keyword_match"] == "true"
        ]
        summaries.append({
            "species": species,
            "canonical_cluster_peaks": len(local_peaks),
            "all_samples": len(species_rows),
            "brain_CNS_keyword_samples": len(brain_rows),
            "developmental_brain_keyword_samples": len(developmental_brain),
            "cluster_detected_all_samples": detected(
                species_rows, "canonical_cluster_detected"
            ),
            "cluster_detected_brain_samples": detected(
                brain_rows, "canonical_cluster_detected"
            ),
            "cluster_detected_developmental_brain_samples": detected(
                developmental_brain, "canonical_cluster_detected"
            ),
            "cluster_detected_nonbrain_samples": detected(
                nonbrain_rows, "canonical_cluster_detected"
            ),
            "brain_classification_method": (
                "transparent sample-description keyword screen"
            ),
        })
    with (OUT / "canonical_fantom5_sample_counts.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(sample_rows[0])
        )
        writer.writeheader()
        writer.writerows(sample_rows)
    summary = {
        "species": summaries,
        "brain_specificity_inference": False,
        "interpretation": (
            "Canonical DNMT3A1 promoter-cluster detection in brain/CNS "
            "samples establishes brain-context activity. Raw detection "
            "fractions do not establish brain enrichment or developmental "
            "specificity."
        ),
    }
    (OUT / "canonical_fantom5_tissue_context_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
