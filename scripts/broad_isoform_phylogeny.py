#!/usr/bin/env python3
"""Summarize phylogenetic distribution of DNMT3A downstream-start products."""

from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/06_isoform_evolution/isoform_species_summary.tsv"
OUT = ROOT / "results/06_isoform_evolution/broad_isoform_phylogeny"
HIGH_ANNOTATION_PRODUCTS = 8


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty output: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def strict_present(row: dict[str, str]) -> bool:
    return int(row["A2_like_products"]) > 0


def broad_present(row: dict[str, str]) -> bool:
    return (
        row["broader_downstream_core_status"]
        != "no_downstream_core_product_detected"
    )


def summarize_group(
    group_type: str, group: str, rows: list[dict[str, str]]
) -> dict[str, object]:
    strict = sum(strict_present(row) for row in rows)
    broad = sum(broad_present(row) for row in rows)
    return {
        "group_type": group_type,
        "group": group,
        "species": len(rows),
        "strict_A2_like_species": strict,
        "strict_A2_like_detection_fraction": strict / len(rows),
        "any_downstream_core_species": broad,
        "any_downstream_core_detection_fraction": broad / len(rows),
        "median_annotated_protein_products": statistics.median(
            int(row["annotated_protein_products"]) for row in rows
        ),
        "interpretation": (
            "annotation-presence distribution; non-detection is not a "
            "biological isoform loss"
        ),
    }


def main() -> None:
    species = read_tsv(INPUT)
    grouped = defaultdict(list)
    for row in species:
        grouped[("mammal_group", row["mammal_group"])].append(row)
        grouped[("order", row["order"])].append(row)
    distribution = [
        summarize_group(group_type, group, rows)
        for (group_type, group), rows in sorted(grouped.items())
    ]

    depth_rows = []
    for label, selected in (
        (
            f"fewer_than_{HIGH_ANNOTATION_PRODUCTS}_products",
            [
                row for row in species
                if int(row["annotated_protein_products"])
                < HIGH_ANNOTATION_PRODUCTS
            ],
        ),
        (
            f"{HIGH_ANNOTATION_PRODUCTS}_or_more_products",
            [
                row for row in species
                if int(row["annotated_protein_products"])
                >= HIGH_ANNOTATION_PRODUCTS
            ],
        ),
    ):
        depth_rows.append(summarize_group("annotation_depth", label, selected))

    evidence_rows = []
    for row in species:
        evidence_rows.append({
            "species": row["species"],
            "order": row["order"],
            "mammal_group": row["mammal_group"],
            "annotated_protein_products": int(
                row["annotated_protein_products"]
            ),
            "strict_A2_like_present": strict_present(row),
            "any_downstream_core_product_present": broad_present(row),
            "broadest_product_status": row[
                "broader_downstream_core_status"
            ],
            "positive_evidence_interpretation": (
                "annotated downstream-core product supports an isoform "
                "hypothesis"
                if broad_present(row) else ""
            ),
            "negative_evidence_interpretation": (
                "" if broad_present(row) else
                "annotation non-detection; biological absence not inferred"
            ),
        })

    eutheria = [row for row in species if row["mammal_group"] == "Eutheria"]
    metatheria = [
        row for row in species if row["mammal_group"] == "Metatheria"
    ]
    prototheria = [
        row for row in species if row["mammal_group"] == "Prototheria"
    ]
    summary = {
        "mammal_species": len(species),
        "orders": len({row["order"] for row in species}),
        "eutherian_orders_with_strict_A2_like_products": len({
            row["order"] for row in eutheria if strict_present(row)
        }),
        "eutherian_orders_total": len({row["order"] for row in eutheria}),
        "eutherian_orders_with_any_downstream_core_product": len({
            row["order"] for row in eutheria if broad_present(row)
        }),
        "metatherian_orders_with_any_downstream_core_product": sorted({
            row["order"] for row in metatheria if broad_present(row)
        }),
        "metatherian_species_with_any_downstream_core_product": sorted(
            row["species"] for row in metatheria if broad_present(row)
        ),
        "prototherian_species_in_package": len(prototheria),
        "prototherian_species_with_any_downstream_core_product": sum(
            broad_present(row) for row in prototheria
        ),
        "annotation_depth_threshold_products": HIGH_ANNOTATION_PRODUCTS,
        "low_depth_broad_detection_fraction": next(
            row["any_downstream_core_detection_fraction"]
            for row in depth_rows
            if row["group"].startswith("fewer")
        ),
        "high_depth_broad_detection_fraction": next(
            row["any_downstream_core_detection_fraction"]
            for row in depth_rows
            if row["group"].startswith(str(HIGH_ANNOTATION_PRODUCTS))
        ),
        "phylogenetic_interpretation": (
            "Downstream-start-compatible products occur broadly across "
            "Eutheria and in two deeply separated marsupial orders. Together "
            "with targeted marsupial transcript evidence, this distribution "
            "is compatible with a therian-ancestral internal-transcript "
            "architecture."
        ),
        "gain_loss_boundary": (
            "The matrix cannot distinguish a therian origin from convergent "
            "eutherian and metatherian recruitment, and it cannot infer "
            "losses from annotation non-detection. Two monotreme "
            "non-detections do not establish absence."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT / "isoform_species_evidence_matrix.tsv", evidence_rows)
    write_tsv(OUT / "isoform_clade_distribution.tsv", distribution)
    write_tsv(OUT / "annotation_depth_sensitivity.tsv", depth_rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
