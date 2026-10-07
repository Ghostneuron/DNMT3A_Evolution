#!/usr/bin/env python3
"""Audit annotation-depth and clade biases in the A2-like product screen."""

from __future__ import annotations

import csv
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/06_isoform_evolution"
SPECIES = RESULTS / "isoform_species_summary.tsv"
PRODUCTS = RESULTS / "isoform_product_classification.tsv"


def depth_bin(count: int) -> str:
    if count <= 2:
        return "1-2_products"
    if count <= 4:
        return "3-4_products"
    if count <= 9:
        return "5-9_products"
    return "10plus_products"


def main() -> None:
    with SPECIES.open(newline="") as handle:
        species_rows = list(csv.DictReader(handle, delimiter="\t"))
    with PRODUCTS.open(newline="") as handle:
        product_rows = list(csv.DictReader(handle, delimiter="\t"))

    strata: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in species_rows:
        strata[("annotation_depth", depth_bin(int(row["annotated_protein_products"])))].append(row)
        strata[("mammal_group", row["mammal_group"])].append(row)
        strata[("order", row["order"])].append(row)

    output = []
    for (stratum_type, stratum), rows in sorted(strata.items()):
        detected = [
            row for row in rows if int(row["A2_like_products"]) > 0
        ]
        output.append({
            "stratum_type": stratum_type,
            "stratum": stratum,
            "species": len(rows),
            "species_with_A2_like_annotation": len(detected),
            "detection_fraction": len(detected) / len(rows),
            "median_annotated_protein_products": statistics.median(
                int(row["annotated_protein_products"]) for row in rows
            ),
            "interpretation": (
                "annotation detection rate; not biological prevalence"
            ),
        })
    fields = list(output[0])
    with (RESULTS / "isoform_annotation_strata.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)

    candidates = [
        row for row in product_rows
        if row["product_class"] == "DNMT3A2_like_downstream_start"
    ]
    leaderless = [
        row for row in product_rows
        if row["product_class"] == "leaderless_downstream_core_product"
    ]
    detected_depth = [
        int(row["annotated_protein_products"]) for row in species_rows
        if int(row["A2_like_products"]) > 0
    ]
    nondetected_depth = [
        int(row["annotated_protein_products"]) for row in species_rows
        if int(row["A2_like_products"]) == 0
    ]
    summary = {
        "species_with_A2_like_annotation": len(detected_depth),
        "species_without_detected_A2_like_annotation": len(nondetected_depth),
        "median_products_detected_species": statistics.median(detected_depth),
        "median_products_nondetected_species": statistics.median(nondetected_depth),
        "candidate_mammal_groups": dict(sorted(Counter(
            row["mammal_group"] for row in candidates
        ).items())),
        "leaderless_downstream_core_mammal_groups": dict(sorted(Counter(
            row["mammal_group"] for row in leaderless
        ).items())),
        "leaderless_non_eutherian_products": [
            {
                "species": row["species"],
                "accession": row["protein_accession"],
                "mammal_group": row["mammal_group"],
                "length_aa": int(row["protein_length_aa"]),
                "full_core_start_0based": int(
                    row["shared_block_full_start_0based"]
                ),
            }
            for row in leaderless if row["mammal_group"] != "Eutheria"
        ],
        "candidate_full_core_start_range_0based": [
            min(int(row["shared_block_full_start_0based"]) for row in candidates),
            max(int(row["shared_block_full_start_0based"]) for row in candidates),
        ],
        "candidate_leader_length_counts": dict(sorted(Counter(
            row["shared_block_product_start_0based"] for row in candidates
        ).items(), key=lambda item: int(item[0]))),
        "gain_loss_inference_ready": False,
        "reason": (
            "Strict A2-like detection is annotation-depth biased; leaderless "
            "downstream-core products require transcript-start validation"
        ),
    }
    (RESULTS / "isoform_annotation_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
