#!/usr/bin/env python3
"""Summarize agreement between independent ESM-2 DNMT3A variant screens."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/07_functional_prioritization/in_silico_esm"
OUTPUT = INPUT / "model_sensitivity"
MODELS = ("esm2_t12_35M_UR50D", "esm2_t30_150M_UR50D")
SCORE = "alternative_vs_reference_log_likelihood_ratio"


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


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average = ((start + 1) + end) / 2
        for index in order[start:end]:
            ranks[index] = average
        start = end
    return ranks


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        raise ValueError("Correlation requires paired vectors of length >= 2")
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum(
        (a - left_mean) * (b - right_mean) for a, b in zip(left, right)
    )
    denominator = math.sqrt(
        sum((a - left_mean) ** 2 for a in left)
        * sum((b - right_mean) ** 2 for b in right)
    )
    return numerator / denominator


def spearman(left: list[float], right: list[float]) -> float:
    return pearson(average_ranks(left), average_ranks(right))


def keyed_scores(path: Path) -> dict[str, dict[str, str]]:
    return {row["variant"]: row for row in read_tsv(path)}


def summarize_set(filename: str, source: str) -> tuple[list[dict[str, object]], dict]:
    model_rows = {
        model: keyed_scores(INPUT / model / filename) for model in MODELS
    }
    shared = sorted(set.intersection(*(set(rows) for rows in model_rows.values())))
    first = model_rows[MODELS[0]]
    second = model_rows[MODELS[1]]
    left = [float(first[variant][SCORE]) for variant in shared]
    right = [float(second[variant][SCORE]) for variant in shared]
    rows: list[dict[str, object]] = []
    for variant, left_score, right_score in zip(shared, left, right):
        left_sign = left_score > 0
        right_sign = right_score > 0
        rows.append({
            "source": source,
            "variant": variant,
            f"{MODELS[0]}_LLR": left_score,
            f"{MODELS[1]}_LLR": right_score,
            "absolute_LLR_difference": abs(left_score - right_score),
            "same_preference_direction": left_sign == right_sign,
            "priority_tier": first[variant].get("priority_tier", ""),
        })
    summary = {
        "source": source,
        "shared_variants": len(shared),
        "pearson_r": pearson(left, right),
        "spearman_rho": spearman(left, right),
        "same_preference_direction_count": sum(
            row["same_preference_direction"] for row in rows
        ),
        "same_preference_direction_fraction": sum(
            row["same_preference_direction"] for row in rows
        ) / len(rows),
        "median_absolute_LLR_difference": sorted(
            float(row["absolute_LLR_difference"]) for row in rows
        )[len(rows) // 2],
    }
    return rows, summary


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    natural_rows, natural_summary = summarize_set(
        "natural_alternative_scores.tsv", "recurrent_natural_alternatives"
    )
    candidate_rows, candidate_summary = summarize_set(
        "candidate_scores.tsv", "candidate_panel"
    )
    write_tsv(OUTPUT / "paired_variant_scores.tsv", natural_rows + candidate_rows)
    tier_one = [
        row for row in candidate_rows if int(row["priority_tier"]) == 1
    ]
    summary = {
        "models": list(MODELS),
        "natural_alternative_concordance": natural_summary,
        "candidate_panel_concordance": candidate_summary,
        "tier_one_candidate_count": len(tier_one),
        "tier_one_same_preference_direction_count": sum(
            row["same_preference_direction"] for row in tier_one
        ),
        "tier_one_variants": [
            {
                "variant": row["variant"],
                "same_preference_direction": row["same_preference_direction"],
                f"{MODELS[0]}_LLR": row[f"{MODELS[0]}_LLR"],
                f"{MODELS[1]}_LLR": row[f"{MODELS[1]}_LLR"],
            }
            for row in tier_one
        ],
        "inference_boundary": (
            "Cross-model agreement shows score robustness to ESM-2 model size; "
            "it does not establish biochemical function, phenotype, or selection."
        ),
    }
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
