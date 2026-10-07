#!/usr/bin/env python3
"""Describe prenatal-to-postnatal shifts in total DNMT3A expression."""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/traits/brain_phenotypes.tsv"
OUTPUT = ROOT / "results/05_brain_integration/expression_phase_contrasts.tsv"


def broad_phase(stage: str) -> str:
    if stage == "embryo":
        return "prenatal"
    if stage:
        return "postnatal_or_later"
    return ""


def main() -> None:
    groups: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    with INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["phenotype"] != "brain_total_DNMT3A_expression":
                continue
            phase = broad_phase(row["developmental_stage"])
            if phase:
                groups[(row["species"], row["brain_region"], phase)].append(
                    math.log2(float(row["value"]) + 1)
                )

    output = []
    keys = sorted({(species, region) for species, region, _ in groups})
    for species, region in keys:
        prenatal = groups.get((species, region, "prenatal"), [])
        postnatal = groups.get((species, region, "postnatal_or_later"), [])
        if not prenatal or not postnatal:
            continue
        prenatal_median = statistics.median(prenatal)
        postnatal_median = statistics.median(postnatal)
        difference = postnatal_median - prenatal_median
        output.append({
            "species": species,
            "brain_region": region,
            "prenatal_samples": len(prenatal),
            "postnatal_or_later_samples": len(postnatal),
            "prenatal_median_log2_CPM_plus_1": prenatal_median,
            "postnatal_median_log2_CPM_plus_1": postnatal_median,
            "postnatal_minus_prenatal_log2": difference,
            "postnatal_to_prenatal_ratio_from_medians": 2 ** difference,
            "direction": "increase" if difference > 0 else "decrease",
            "interpretation_limit": (
                "Descriptive within-species contrast; age distribution and "
                "anatomical definition differ across phases"
            ),
        })

    fields = [
        "species", "brain_region", "prenatal_samples",
        "postnatal_or_later_samples", "prenatal_median_log2_CPM_plus_1",
        "postnatal_median_log2_CPM_plus_1",
        "postnatal_minus_prenatal_log2",
        "postnatal_to_prenatal_ratio_from_medians", "direction",
        "interpretation_limit",
    ]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    decreases = sum(row["direction"] == "decrease" for row in output)
    print(
        f"Wrote {len(output)} contrasts to {OUTPUT}; "
        f"{decreases}/{len(output)} are decreases"
    )


if __name__ == "__main__":
    main()
