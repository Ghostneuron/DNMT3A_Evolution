#!/usr/bin/env python3
"""Analyze independent MeCP2 gene sets against DNMT3A isoform mCA effects."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact


SET_DEFINITIONS = {
    "core_MR": (lambda row: row["Core"] == "MR", lambda row: row["Core"] == "Unchanged"),
    "core_MA": (lambda row: row["Core"] == "MA", lambda row: row["Core"] == "Unchanged"),
    "L4_MR": (lambda row: row["L4"] == "MR", lambda row: row["L4"] == "Unchanged"),
    "L5_MR": (lambda row: row["L5"] == "MR", lambda row: row["L5"] == "Unchanged"),
    "PV_MR": (lambda row: row["PV"] == "MR", lambda row: row["PV"] == "Unchanged"),
    "SST_MR": (lambda row: row["SST"] == "MR", lambda row: row["SST"] == "Unchanged"),
    "any_subclass_MR": (
        lambda row: int(row["MR_subclass_count"]) >= 1,
        lambda row: row["all_subclasses_unchanged"] == "True",
    ),
    "recurrent_subclass_MR": (
        lambda row: int(row["MR_subclass_count"]) >= 2,
        lambda row: row["all_subclasses_unchanged"] == "True",
    ),
}


def bh_adjust(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [1.0] * len(values)
    running = 1.0
    for rank_index in range(len(order) - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[index] * len(values) / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def interval(values: list[float]) -> tuple[float, float]:
    low, high = np.percentile(values, [2.5, 97.5])
    return float(low), float(high)


def bootstrap_median_difference(
    member: list[dict[str, object]],
    control: list[dict[str, object]],
    field: str,
    bootstraps: int,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    observed = statistics.median(float(row[field]) for row in member) - statistics.median(
        float(row[field]) for row in control
    )
    chromosomes = sorted({str(row["chrom"]) for row in member + control})
    member_by_chrom = {
        chrom: np.array([float(row[field]) for row in member if row["chrom"] == chrom])
        for chrom in chromosomes
    }
    control_by_chrom = {
        chrom: np.array([float(row[field]) for row in control if row["chrom"] == chrom])
        for chrom in chromosomes
    }
    estimates = []
    for _ in range(bootstraps):
        selected = rng.choice(chromosomes, len(chromosomes), replace=True)
        sampled_member = np.concatenate(
            [member_by_chrom[chrom] for chrom in selected if len(member_by_chrom[chrom])]
        )
        sampled_control = np.concatenate(
            [control_by_chrom[chrom] for chrom in selected if len(control_by_chrom[chrom])]
        )
        if len(sampled_member) and len(sampled_control):
            estimates.append(float(np.median(sampled_member) - np.median(sampled_control)))
    low, high = interval(estimates)
    return observed, low, high


def build(
    methylation_path: Path,
    expression_path: Path,
    gene_class_path: Path,
    output_dir: Path,
    bootstraps: int,
    seed: int,
) -> dict[str, object]:
    with gene_class_path.open(newline="") as handle:
        classes = {row["gene"]: row for row in csv.DictReader(handle, delimiter="\t")}
    expression = {}
    with gzip.open(expression_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expression[(row["contrast"], row["gene"])] = row
    grouped = {}
    with gzip.open(methylation_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene_class = classes.get(row["gene"])
            expr = expression.get((row["contrast"], row["gene"]))
            if gene_class is None or expr is None:
                continue
            grouped.setdefault(row["contrast"], []).append({**row, **expr, **gene_class})

    rng = np.random.default_rng(seed)
    results = []
    pvalues = []
    for contrast in sorted(grouped):
        rows = grouped[contrast]
        threshold = float(np.percentile([float(row["absolute_difference"]) for row in rows], 10))
        for name, (is_member, is_control) in SET_DEFINITIONS.items():
            member = [row for row in rows if is_member(row)]
            control = [row for row in rows if is_control(row)]
            if not member or not control:
                continue
            metrics = {}
            for field, label, eligible in (
                ("wt_mean_corrected_mCA", "WT_mCA", lambda row: True),
                ("absolute_difference", "mCA_effect", lambda row: True),
                ("relative_to_WT", "mCA_retention", lambda row: True),
                (
                    "relative_to_WT",
                    "mCA_retention_WTge0p005",
                    lambda row: float(row["wt_mean_corrected_mCA"]) >= 0.005,
                ),
                ("expression_difference", "expression_effect", lambda row: True),
            ):
                metric_member = [
                    row for row in member if row[field] != "" and eligible(row)
                ]
                metric_control = [
                    row for row in control if row[field] != "" and eligible(row)
                ]
                observed, low, high = bootstrap_median_difference(
                    metric_member, metric_control, field, bootstraps, rng
                )
                metrics[f"{label}_median_difference"] = observed
                metrics[f"{label}_block_interval_low"] = low
                metrics[f"{label}_block_interval_high"] = high
                metrics[f"{label}_member_genes"] = len(metric_member)
                metrics[f"{label}_control_genes"] = len(metric_control)
            member_bottom = sum(float(row["absolute_difference"]) <= threshold for row in member)
            control_bottom = sum(float(row["absolute_difference"]) <= threshold for row in control)
            odds, pvalue = fisher_exact(
                [
                    [member_bottom, len(member) - member_bottom],
                    [control_bottom, len(control) - control_bottom],
                ]
            )
            results.append(
                {
                    "contrast": contrast,
                    "gene_set": name,
                    "member_genes": len(member),
                    "control_genes": len(control),
                    **metrics,
                    "most_negative_mCA_decile_threshold": threshold,
                    "member_in_most_negative_decile": member_bottom,
                    "control_in_most_negative_decile": control_bottom,
                    "most_negative_decile_odds_ratio": float(odds),
                    "fisher_p": float(pvalue),
                }
            )
            pvalues.append(float(pvalue))
    for row, qvalue in zip(results, bh_adjust(pvalues)):
        row["fisher_BH_q"] = qvalue

    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "mecp2_gene_set_enrichment.tsv"
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    report = {
        "complete": bool(results),
        "contrasts": len(grouped),
        "tests": len(results),
        "chromosome_block_bootstraps": bootstraps,
        "output": str(table_path),
        "interpretation_boundary": (
            "Gene-set association with chromosome-block robustness; not animal-level "
            "inference and not evidence that mCA change caused expression change."
        ),
    }
    (output_dir / "mecp2_gene_set_enrichment.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--methylation", type=Path,
        default=Path("results/mch_gene_body_comparison/gene_body_mCA_contrasts.tsv.gz"),
    )
    parser.add_argument(
        "--expression", type=Path,
        default=Path("results/mch_gene_body_comparison/cortex_expression_contrasts.tsv.gz"),
    )
    parser.add_argument(
        "--gene-classes", type=Path,
        default=Path("results/mch_gene_body_comparison/moore_2025_mecp2_gene_classes.tsv"),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/mch_gene_body_comparison")
    )
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=161803)
    args = parser.parse_args()
    print(json.dumps(build(args.methylation, args.expression, args.gene_classes, args.output_dir, args.bootstraps, args.seed), indent=2))


if __name__ == "__main__":
    main()
