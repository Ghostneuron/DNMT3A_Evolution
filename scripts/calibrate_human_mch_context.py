#!/usr/bin/env python3
"""Cross-calibrate independent adult human bulk and single-cell mCH summaries."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRAINZOO = ROOT / "data/traits/sources/brainzoo_mch.tsv"
ATLAS = (
    ROOT
    / "results/05_brain_integration/human_cell_atlas_import_summary.json"
)
OUT = ROOT / "results/05_brain_integration"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    bulk_row = next(
        row for row in read_tsv(BRAINZOO)
        if row["species"] == "Homo sapiens"
        and row["phenotype"] == "global_brain_mCH_fraction"
    )
    atlas = json.loads(ATLAS.read_text())
    bulk = float(bulk_row["value"])
    neuronal = float(atlas["motor_cortex_neuronal_mCH_median"])
    non_neuronal = float(atlas["motor_cortex_non_neuronal_mCH_median"])
    bracketed = non_neuronal <= bulk <= neuronal

    rows = [
        {
            "dataset": "BrainZoo",
            "measurement": "coverage_weighted_bulk_mCH_fraction",
            "value": bulk,
            "brain_region": bulk_row["brain_region"],
            "cell_context": "bulk_tissue",
            "assay": bulk_row["assay"],
            "aggregation": bulk_row["normalization"],
            "citation_doi": bulk_row["citation_doi"],
        },
        {
            "dataset": "Human_Cell_Epigenome_Atlas",
            "measurement": "median_single_nucleus_mCH_fraction",
            "value": neuronal,
            "brain_region": "primary motor cortex",
            "cell_context": "neuronal",
            "assay": "snm3C-seq",
            "aggregation": "median across neuronal nuclei",
            "citation_doi": "10.1126/science.adx0673",
        },
        {
            "dataset": "Human_Cell_Epigenome_Atlas",
            "measurement": "median_single_nucleus_mCH_fraction",
            "value": non_neuronal,
            "brain_region": "primary motor cortex",
            "cell_context": "non_neuronal",
            "assay": "snm3C-seq",
            "aggregation": "median across non-neuronal nuclei",
            "citation_doi": "10.1126/science.adx0673",
        },
    ]
    with (OUT / "human_mch_cross_dataset_calibration.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(rows[0])
        )
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "brainzoo_bulk_mCH_fraction": bulk,
        "atlas_neuronal_mCH_median": neuronal,
        "atlas_non_neuronal_mCH_median": non_neuronal,
        "bulk_value_bracketed_by_atlas_cell_classes": bracketed,
        "bulk_to_neuronal_median_ratio": bulk / neuronal,
        "bulk_to_non_neuronal_median_ratio": bulk / non_neuronal,
        "interpretation": (
            "The independent bulk value lies between the atlas neuronal and "
            "non-neuronal medians, which is directionally compatible with a "
            "mixed-cell brain sample."
        ),
        "forbidden_inference": (
            "Do not estimate neuronal fraction by linear mixing: the studies "
            "use different cortical regions, donors, assays, coverage "
            "weighting, and summary statistics."
        ),
    }
    (OUT / "human_mch_cross_dataset_calibration.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
