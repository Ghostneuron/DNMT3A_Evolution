#!/usr/bin/env python3
"""Summarize developmental total-DNMT3A expression without pooling tissues."""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/traits/brain_phenotypes.tsv"
OUTPUT = ROOT / "results/05_brain_integration/expression_stage_summary.tsv"


def main() -> None:
    with INPUT.open(newline="") as handle:
        rows = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["phenotype"] == "brain_total_DNMT3A_expression"
        ]
    groups: dict[tuple[str, str, str, str, str, str], list[float]] = defaultdict(list)
    for row in rows:
        key = (
            row["species"], row["brain_region"],
            row["developmental_stage"], row["age_value"], row["age_unit"],
            row["unit"],
        )
        groups[key].append(float(row["value"]))

    output = []
    for (species, region, stage, age_value, age_unit, unit), values in groups.items():
        output.append({
            "species": species,
            "brain_region": region,
            "developmental_stage": stage,
            "age_value": age_value,
            "age_unit": age_unit,
            "unit": unit,
            "replicates": len(values),
            "mean": statistics.fmean(values),
            "median": statistics.median(values),
            "minimum": min(values),
            "maximum": max(values),
        })
    output.sort(key=lambda row: (
        row["species"], row["brain_region"], row["developmental_stage"],
        row["age_unit"], row["age_value"],
    ))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "species", "brain_region", "developmental_stage", "age_value",
        "age_unit", "unit",
        "replicates", "mean", "median", "minimum", "maximum",
    ]
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    print(f"Wrote {len(output)} tissue-stage summaries to {OUTPUT}")


if __name__ == "__main__":
    main()
