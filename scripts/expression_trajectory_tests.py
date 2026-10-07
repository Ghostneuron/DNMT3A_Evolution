#!/usr/bin/env python3
"""Permutation-test monotonic total-DNMT3A developmental trajectories."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/traits/brain_phenotypes.tsv"
RESULTS = ROOT / "results/05_brain_integration"
OUTPUT = RESULTS / "expression_trajectory_monotonicity.tsv"
PERMUTATIONS = 20_000

POSTNATAL_ORDER = {
    "newborn": 0,
    "infant": 1,
    "toddler": 2,
    "school": 3,
    "teenager": 4,
    "youngAdult": 5,
    "youngMidAge": 6,
    "olderMidAge": 7,
    "senior": 8,
}


def rankdata(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    position = 0
    while position < len(values):
        end = position + 1
        while end < len(values) and values[order[end]] == values[order[position]]:
            end += 1
        average_rank = (position + 1 + end) / 2
        for index in order[position:end]:
            ranks[index] = average_rank
        position = end
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    x_mean = statistics.fmean(x)
    y_mean = statistics.fmean(y)
    numerator = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y))
    x_ss = sum((a - x_mean) ** 2 for a in x)
    y_ss = sum((b - y_mean) ** 2 for b in y)
    if x_ss == 0 or y_ss == 0:
        return 0.0
    return numerator / math.sqrt(x_ss * y_ss)


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(rankdata(x), rankdata(y))


def permutation_p(
    x: list[float], y: list[float], permutations: int, seed: int
) -> float:
    observed = abs(spearman(x, y))
    rng = random.Random(seed)
    shuffled = y.copy()
    extreme = 0
    for _ in range(permutations):
        rng.shuffle(shuffled)
        if abs(spearman(x, shuffled)) >= observed - 1e-12:
            extreme += 1
    return (extreme + 1) / (permutations + 1)


def age_scale(value: str, unit: str) -> float:
    number = float(value)
    unit = unit.lower()
    if unit.startswith("day"):
        return number
    if unit.startswith("week"):
        return number * 7
    if unit.startswith("month"):
        return number * 30.4375
    if unit.startswith("year"):
        return number * 365.25
    return number


def stage_sort_key(row: dict[str, str]) -> tuple[int, float]:
    stage = row["developmental_stage"]
    phase = 0 if stage == "embryo" else 1
    if row["age_value"]:
        return phase, age_scale(row["age_value"], row["age_unit"])
    return phase, float(POSTNATAL_ORDER.get(stage, 100))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    count = len(p_values)
    order = sorted(range(count), key=p_values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for reverse_rank, index in enumerate(reversed(order), start=1):
        rank = count - reverse_rank + 1
        running = min(running, p_values[index] * count / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def main() -> None:
    samples: dict[tuple[str, str, tuple[str, str, str]], list[float]] = defaultdict(list)
    stage_rows: dict[tuple[str, str, tuple[str, str, str]], dict[str, str]] = {}
    with INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["phenotype"] != "brain_total_DNMT3A_expression":
                continue
            stage_id = (
                row["developmental_stage"], row["age_value"], row["age_unit"]
            )
            key = (row["species"], row["brain_region"], stage_id)
            samples[key].append(math.log2(float(row["value"]) + 1))
            stage_rows[key] = row

    grouped: dict[tuple[str, str], list[tuple[dict[str, str], list[float]]]] = (
        defaultdict(list)
    )
    for (species, region, _), values in samples.items():
        grouped[(species, region)].append(
            (stage_rows[(species, region, _)], values)
        )

    output: list[dict[str, object]] = []
    for (species, region), stages in sorted(grouped.items()):
        stages.sort(key=lambda item: stage_sort_key(item[0]))
        x = list(range(len(stages)))
        y = [statistics.median(values) for _, values in stages]
        rho = spearman(x, y)
        seed = int.from_bytes(
            hashlib.sha256(f"{species}|{region}".encode()).digest()[:8],
            "big",
        )
        p_value = permutation_p(x, y, PERMUTATIONS, seed)
        output.append({
            "species": species,
            "brain_region": region,
            "ordered_stages": len(stages),
            "samples": sum(len(values) for _, values in stages),
            "spearman_rho_stage_vs_median_log2_CPM": rho,
            "permutation_p_two_sided": p_value,
            "permutations": PERMUTATIONS,
            "direction": "decrease" if rho < 0 else "increase",
            "stage_scale": "within-species ordinal developmental stage",
            "interpretation_limit": (
                "Descriptive gene-level trajectory; stage spacing and anatomy "
                "are not identical across species"
            ),
        })

    adjusted = benjamini_hochberg(
        [float(row["permutation_p_two_sided"]) for row in output]
    )
    for row, q_value in zip(output, adjusted):
        row["BH_q_across_12_trajectories"] = q_value

    fields = [
        "species", "brain_region", "ordered_stages", "samples",
        "spearman_rho_stage_vs_median_log2_CPM",
        "permutation_p_two_sided", "BH_q_across_12_trajectories",
        "permutations", "direction", "stage_scale", "interpretation_limit",
    ]
    RESULTS.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    summary = {
        "trajectories": len(output),
        "permutations_per_trajectory": PERMUTATIONS,
        "negative_rho": sum(float(row["spearman_rho_stage_vs_median_log2_CPM"]) < 0 for row in output),
        "BH_q_below_0.05": sum(float(row["BH_q_across_12_trajectories"]) < 0.05 for row in output),
    }
    (RESULTS / "expression_trajectory_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
