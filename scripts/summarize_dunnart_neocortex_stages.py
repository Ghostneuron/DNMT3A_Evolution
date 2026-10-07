#!/usr/bin/env python3
"""Summarize DNMT3A junction evidence across dunnart neocortex stages."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/06_isoform_evolution"
PATTERN = "dunnart_P*_neocortex_rep*_junction_summary.json"
COMPONENTS = (
    "full_length_first_junction",
    "internal_067_first_junction",
    "internal_068_first_junction",
    "common_downstream_junction",
)


def stage_from_sample(sample: str) -> str:
    match = re.search(r"\b(P\d+)\b", sample)
    if not match:
        raise ValueError(f"Cannot determine stage from sample: {sample}")
    return match.group(1)


def exact_permutation_p_greater(
    first: list[float],
    second: list[float],
) -> float | None:
    """Exact animal-label permutation P for mean(first) > mean(second)."""
    if not first or len(first) != len(second):
        return None
    pooled = first + second
    group_size = len(first)
    observed = sum(first) / group_size - sum(second) / group_size
    exceedances = 0
    permutations = 0
    indices = set(range(len(pooled)))
    for selected_tuple in combinations(range(len(pooled)), group_size):
        selected = set(selected_tuple)
        complement = indices - selected
        statistic = (
            sum(pooled[index] for index in selected) / group_size
            - sum(pooled[index] for index in complement) / group_size
        )
        exceedances += statistic >= observed - 1e-15
        permutations += 1
    return exceedances / permutations


def summarize_stages(records: list[dict]) -> tuple[list[dict], dict]:
    per_stage_counts: dict[str, Counter] = defaultdict(Counter)
    per_stage_pairs = Counter()
    per_stage_replicates = Counter()
    per_stage_internal_rates = defaultdict(list)
    per_stage_internal_common_ratios = defaultdict(list)
    per_stage_internal_full_ratios = defaultdict(list)
    replicate_rows = []
    for record in records:
        stage = stage_from_sample(record["sample"])
        counts = record["junction_supporting_pair_counts"]
        pair_count = record["paired_read_count"]
        row = {
            "stage": stage,
            "run_accession": record["run_accession"],
            "sample": record["sample"],
            "paired_read_count": pair_count,
        }
        for component in COMPONENTS:
            count = counts[component]
            row[f"{component}_pairs"] = count
            row[f"{component}_pairs_per_million"] = round(
                count * 1_000_000 / pair_count, 4
            )
            per_stage_counts[stage][component] += count
        diagnostic = (
            counts["full_length_first_junction"]
            + counts["internal_067_first_junction"]
            + counts["internal_068_first_junction"]
        )
        row["internal_067_fraction_of_diagnostic_first_junctions"] = (
            round(counts["internal_067_first_junction"] / diagnostic, 4)
            if diagnostic else None
        )
        row["internal_067_to_common_downstream_pair_ratio"] = (
            round(
                counts["internal_067_first_junction"]
                / counts["common_downstream_junction"],
                4,
            )
            if counts["common_downstream_junction"] else None
        )
        per_stage_internal_rates[stage].append(
            counts["internal_067_first_junction"] * 1_000_000 / pair_count
        )
        if counts["common_downstream_junction"]:
            per_stage_internal_common_ratios[stage].append(
                counts["internal_067_first_junction"]
                / counts["common_downstream_junction"]
            )
        if counts["full_length_first_junction"]:
            per_stage_internal_full_ratios[stage].append(
                counts["internal_067_first_junction"]
                / counts["full_length_first_junction"]
            )
        replicate_rows.append(row)
        per_stage_pairs[stage] += pair_count
        per_stage_replicates[stage] += 1

    stage_summary = {}
    for stage in sorted(per_stage_pairs, key=lambda value: int(value[1:])):
        pair_count = per_stage_pairs[stage]
        counts = per_stage_counts[stage]
        diagnostic = sum(
            counts[component]
            for component in (
                "full_length_first_junction",
                "internal_067_first_junction",
                "internal_068_first_junction",
            )
        )
        stage_summary[stage] = {
            "replicate_count": per_stage_replicates[stage],
            "pooled_paired_read_count": pair_count,
            "junction_supporting_pair_counts": {
                component: counts[component] for component in COMPONENTS
            },
            "pairs_per_million": {
                component: round(
                    counts[component] * 1_000_000 / pair_count, 4
                )
                for component in COMPONENTS
            },
            "internal_067_fraction_of_diagnostic_first_junctions": (
                round(counts["internal_067_first_junction"] / diagnostic, 4)
                if diagnostic else None
            ),
            "internal_067_to_common_downstream_pair_ratio": (
                round(
                    counts["internal_067_first_junction"]
                    / counts["common_downstream_junction"],
                    4,
                )
                if counts["common_downstream_junction"] else None
            ),
            "internal_067_to_full_length_first_pair_ratio": (
                round(
                    counts["internal_067_first_junction"]
                    / counts["full_length_first_junction"],
                    4,
                )
                if counts["full_length_first_junction"] else None
            ),
            "replicate_internal_pairs_per_million_range": [
                round(min(per_stage_internal_rates[stage]), 4),
                round(max(per_stage_internal_rates[stage]), 4),
            ],
            "replicate_internal_to_common_ratio_range": [
                round(min(per_stage_internal_common_ratios[stage]), 4),
                round(max(per_stage_internal_common_ratios[stage]), 4),
            ],
            "replicate_internal_to_full_length_ratio_range": [
                round(min(per_stage_internal_full_ratios[stage]), 4),
                round(max(per_stage_internal_full_ratios[stage]), 4),
            ],
        }
    contrasts = {}
    if "P12" in stage_summary and "P20" in stage_summary:
        p12 = stage_summary["P12"]
        p20 = stage_summary["P20"]
        contrasts["P12_over_P20"] = {
            "internal_pairs_per_million_fold": round(
                p12["pairs_per_million"]["internal_067_first_junction"]
                / p20["pairs_per_million"]["internal_067_first_junction"],
                4,
            ),
            "full_length_pairs_per_million_fold": round(
                p12["pairs_per_million"]["full_length_first_junction"]
                / p20["pairs_per_million"]["full_length_first_junction"],
                4,
            ),
            "common_downstream_pairs_per_million_fold": round(
                p12["pairs_per_million"]["common_downstream_junction"]
                / p20["pairs_per_million"]["common_downstream_junction"],
                4,
            ),
            "internal_to_common_ratio_fold": round(
                p12["internal_067_to_common_downstream_pair_ratio"]
                / p20["internal_067_to_common_downstream_pair_ratio"],
                4,
            ),
            "internal_to_full_length_ratio_fold": round(
                p12["internal_067_to_full_length_first_pair_ratio"]
                / p20["internal_067_to_full_length_first_pair_ratio"],
                4,
            ),
            "all_P12_internal_rates_exceed_all_P20_rates": (
                min(per_stage_internal_rates["P12"])
                > max(per_stage_internal_rates["P20"])
            ),
            "all_P12_internal_to_common_ratios_exceed_all_P20_ratios": (
                min(per_stage_internal_common_ratios["P12"])
                > max(per_stage_internal_common_ratios["P20"])
            ),
            "all_P12_internal_to_full_length_ratios_exceed_all_P20_ratios": (
                min(per_stage_internal_full_ratios["P12"])
                > max(per_stage_internal_full_ratios["P20"])
            ),
            "animal_level_exact_permutation_p_greater": {
                "internal_pairs_per_million": exact_permutation_p_greater(
                    per_stage_internal_rates["P12"],
                    per_stage_internal_rates["P20"],
                ),
                "internal_to_common_ratio": exact_permutation_p_greater(
                    per_stage_internal_common_ratios["P12"],
                    per_stage_internal_common_ratios["P20"],
                ),
                "internal_to_full_length_ratio": exact_permutation_p_greater(
                    per_stage_internal_full_ratios["P12"],
                    per_stage_internal_full_ratios["P20"],
                ),
            },
        }
    replicated_stages = (
        {"P12", "P20"}.issubset(per_stage_replicates)
        and all(per_stage_replicates[stage] >= 2 for stage in ("P12", "P20"))
    )
    summary = {
        "stages": stage_summary,
        "descriptive_contrasts": contrasts,
        "stage_replication_sufficient_for_directional_claim": replicated_stages,
        "interpretation_rule": (
            "With at least two animals per stage, a replicated directional "
            "pattern may be stated. Exact label-permutation tests treat "
            "animals, not read pairs, as replicates and remain low-powered "
            "with two or three animals per stage."
            if replicated_stages else
            "Stage comparisons are descriptive until each stage has at least "
            "two biological replicates; read pairs are not biological replicates."
        ),
    }
    return replicate_rows, summary


def main() -> None:
    paths = sorted(RESULTS.glob(PATTERN))
    if not paths:
        raise SystemExit(f"No inputs matching {RESULTS / PATTERN}")
    records = [json.loads(path.read_text()) for path in paths]
    rows, summary = summarize_stages(records)
    output = RESULTS / "dunnart_neocortex_stage_junctions.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (RESULTS / "dunnart_neocortex_stage_junctions_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
