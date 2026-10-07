#!/usr/bin/env python3
"""Relate cortex chromatin context to DNMT3A isoform-dependent gene-body mCA."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np

try:
    from scripts.model_mch_gene_body_covariates import EXPRESSION_COLUMN, ols, zscore
except ModuleNotFoundError:  # Direct execution adds scripts/, not the project root.
    from model_mch_gene_body_covariates import EXPRESSION_COLUMN, ols, zscore


FEATURES = (
    "dnmt3a_flag_gene_body_enrichment",
    "dnmt3a_flag_promoter_enrichment",
    "h3k4me3_gene_body_enrichment",
    "h3k4me3_promoter_enrichment",
    "h3k27me3_gene_body_enrichment",
    "h3k27me3_promoter_enrichment",
)
JOINT_FEATURES = (
    "dnmt3a_flag_gene_body_enrichment",
    "h3k27me3_gene_body_enrichment",
)


def add_enrichment_features(row: dict[str, str]) -> dict[str, float]:
    values: dict[str, float] = {}
    for track in ("dnmt3a_flag", "h3k4me3", "h3k27me3"):
        for region in ("gene_body", "promoter"):
            signal = float(row[f"{track}_{region}_mean"])
            control = float(row[f"input_{region}_mean"])
            values[f"{track}_{region}_enrichment"] = float(
                np.log1p(signal) - np.log1p(control)
            )
    return values


def analyze_feature(
    rows: list[dict[str, object]], feature: str, bootstraps: int, seed: int
) -> dict[str, object]:
    y = zscore(np.array([float(row["absolute_difference"]) for row in rows]))
    lengths = np.array([float(row["length_bp"]) for row in rows])
    base = [
        zscore(np.array([float(row["wt_mean_corrected_mCA"]) for row in rows])),
        zscore(np.log10(lengths)),
        zscore(
            np.log10(
                np.array([float(row["minimum_sites_across_samples"]) for row in rows])
                / lengths
            )
        ),
        zscore(np.array([float(row["matched_WT_log2_CPM"]) for row in rows])),
    ]
    focal = zscore(np.array([float(row[feature]) for row in rows]))
    chrom = np.array([str(row["chrom"]) for row in rows])
    _, reduced_r2 = ols(y, base)
    extended, extended_r2 = ols(y, [*base, focal])
    partial_r2 = max(0.0, (extended_r2 - reduced_r2) / max(1e-15, 1.0 - reduced_r2))

    blocks = sorted(set(chrom))
    rng = np.random.default_rng(seed)
    extended_design = np.column_stack([np.ones(len(y)), *base, focal])
    block_stats = {}
    for block in blocks:
        indices = np.flatnonzero(chrom == block)
        design = extended_design[indices]
        outcome = y[indices]
        block_stats[block] = (
            design.T @ design,
            design.T @ outcome,
            float(outcome @ outcome),
            float(outcome.sum()),
            len(outcome),
        )

    def fit_from_stats(selected: np.ndarray, columns: int) -> tuple[np.ndarray, float]:
        xtx = sum((block_stats[str(block)][0][:columns, :columns] for block in selected))
        xty = sum((block_stats[str(block)][1][:columns] for block in selected))
        yty = sum(block_stats[str(block)][2] for block in selected)
        ysum = sum(block_stats[str(block)][3] for block in selected)
        count = sum(block_stats[str(block)][4] for block in selected)
        beta = np.linalg.lstsq(xtx, xty, rcond=None)[0]
        sse = max(0.0, float(yty - beta @ xty))
        total = max(1e-15, float(yty - ysum * ysum / count))
        return beta, 1.0 - sse / total

    coefficients: list[float] = []
    partials: list[float] = []
    for _ in range(bootstraps):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        _, reduced_boot_r2 = fit_from_stats(selected, 5)
        extended_boot, extended_boot_r2 = fit_from_stats(selected, 6)
        coefficients.append(float(extended_boot[-1]))
        partials.append(
            max(
                0.0,
                (extended_boot_r2 - reduced_boot_r2)
                / max(1e-15, 1.0 - reduced_boot_r2),
            )
        )
    coefficient_low, coefficient_high = np.percentile(coefficients, [2.5, 97.5])
    partial_low, partial_high = np.percentile(partials, [2.5, 97.5])
    return {
        "genes": len(rows),
        "chromosome_blocks": len(blocks),
        "standardized_feature_beta": float(extended[-1]),
        "feature_beta_block_interval_low": float(coefficient_low),
        "feature_beta_block_interval_high": float(coefficient_high),
        "base_model_r_squared": reduced_r2,
        "extended_model_r_squared": extended_r2,
        "feature_partial_r_squared": partial_r2,
        "partial_r_squared_block_interval_low": float(partial_low),
        "partial_r_squared_block_interval_high": float(partial_high),
        "chromosome_block_bootstraps": bootstraps,
    }


def analyze_joint_features(
    rows: list[dict[str, object]], features: tuple[str, str], bootstraps: int, seed: int
) -> dict[str, object]:
    y = zscore(np.array([float(row["absolute_difference"]) for row in rows]))
    lengths = np.array([float(row["length_bp"]) for row in rows])
    base = [
        zscore(np.array([float(row["wt_mean_corrected_mCA"]) for row in rows])),
        zscore(np.log10(lengths)),
        zscore(np.log10(np.array([float(row["minimum_sites_across_samples"]) for row in rows]) / lengths)),
        zscore(np.array([float(row["matched_WT_log2_CPM"]) for row in rows])),
    ]
    focal = [zscore(np.array([float(row[feature]) for row in rows])) for feature in features]
    chrom = np.array([str(row["chrom"]) for row in rows])
    _, reduced_r2 = ols(y, base)
    extended, extended_r2 = ols(y, [*base, *focal])
    blocks = sorted(set(chrom))
    design = np.column_stack([np.ones(len(y)), *base, *focal])
    stats = {}
    for block in blocks:
        indices = np.flatnonzero(chrom == block)
        x, outcome = design[indices], y[indices]
        stats[block] = (x.T @ x, x.T @ outcome, float(outcome @ outcome), float(outcome.sum()), len(outcome))
    rng = np.random.default_rng(seed)
    coefficients = [[], []]
    partials = []
    for _ in range(bootstraps):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        xtx = sum(stats[str(block)][0] for block in selected)
        xty = sum(stats[str(block)][1] for block in selected)
        yty = sum(stats[str(block)][2] for block in selected)
        ysum = sum(stats[str(block)][3] for block in selected)
        count = sum(stats[str(block)][4] for block in selected)
        beta = np.linalg.lstsq(xtx, xty, rcond=None)[0]
        total = max(1e-15, yty - ysum * ysum / count)
        full_r2 = 1.0 - max(0.0, yty - float(beta @ xty)) / total
        base_beta = np.linalg.lstsq(xtx[:5, :5], xty[:5], rcond=None)[0]
        base_r2 = 1.0 - max(0.0, yty - float(base_beta @ xty[:5])) / total
        coefficients[0].append(float(beta[-2]))
        coefficients[1].append(float(beta[-1]))
        partials.append(max(0.0, (full_r2 - base_r2) / max(1e-15, 1.0 - base_r2)))
    intervals = [np.percentile(values, [2.5, 97.5]) for values in coefficients]
    partial_interval = np.percentile(partials, [2.5, 97.5])
    return {
        "genes": len(rows),
        "chromosome_blocks": len(blocks),
        "dnmt3a_flag_standardized_beta": float(extended[-2]),
        "dnmt3a_flag_beta_block_interval_low": float(intervals[0][0]),
        "dnmt3a_flag_beta_block_interval_high": float(intervals[0][1]),
        "h3k27me3_standardized_beta": float(extended[-1]),
        "h3k27me3_beta_block_interval_low": float(intervals[1][0]),
        "h3k27me3_beta_block_interval_high": float(intervals[1][1]),
        "base_model_r_squared": reduced_r2,
        "joint_model_r_squared": extended_r2,
        "joint_partial_r_squared": max(0.0, (extended_r2 - reduced_r2) / max(1e-15, 1.0 - reduced_r2)),
        "joint_partial_r_squared_block_interval_low": float(partial_interval[0]),
        "joint_partial_r_squared_block_interval_high": float(partial_interval[1]),
        "chromosome_block_bootstraps": bootstraps,
    }


def build(
    methylation_path: Path,
    expression_path: Path,
    chromatin_path: Path,
    output_dir: Path,
    bootstraps: int,
    seed: int,
) -> dict[str, object]:
    with expression_path.open(newline="") as handle:
        expression = {row["gene"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with gzip.open(chromatin_path, "rt", newline="") as handle:
        chromatin = {
            (row["chrom"], row["start"], row["end"], row["gene"]): add_enrichment_features(row)
            for row in csv.DictReader(handle, delimiter="\t")
        }

    grouped: dict[str, list[dict[str, object]]] = {}
    with gzip.open(methylation_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expr = expression.get(row["gene"])
            features = chromatin.get((row["chrom"], row["start"], row["end"], row["gene"]))
            if expr is None or features is None:
                continue
            row["matched_WT_log2_CPM"] = expr[EXPRESSION_COLUMN[row["contrast"]]]
            grouped.setdefault(row["contrast"], []).append({**row, **features})

    results = []
    joint_results = []
    iteration = 0
    for contrast in sorted(grouped):
        all_rows = grouped[contrast]
        for sensitivity, rows in (
            ("all_genes", all_rows),
            (
                "WT_mCA_ge_0.005",
                [row for row in all_rows if float(row["wt_mean_corrected_mCA"]) >= 0.005],
            ),
        ):
            for feature in FEATURES:
                result = analyze_feature(rows, feature, bootstraps, seed + iteration)
                results.append(
                    {"contrast": contrast, "sensitivity": sensitivity, "feature": feature, **result}
                )
                iteration += 1
            joint_results.append(
                {
                    "contrast": contrast,
                    "sensitivity": sensitivity,
                    **analyze_joint_features(rows, JOINT_FEATURES, bootstraps, seed + iteration),
                }
            )
            iteration += 1

    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "mch_chromatin_context_models.tsv"
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    joint_table_path = output_dir / "mch_dnmt3a_h3k27me3_joint_models.tsv"
    with joint_table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(joint_results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(joint_results)
    report = {
        "complete": bool(results),
        "output": str(table_path),
        "contrasts": len(grouped),
        "models": len(results),
        "joint_models": len(joint_results),
        "joint_output": str(joint_table_path),
        "feature_definition": "log1p(ChIP mean signal) minus log1p(input mean signal)",
        "outcome": "standardized genotype-minus-WT corrected gene-body mCA",
        "base_predictors": [
            "WT corrected gene-body mCA",
            "log10 gene length",
            "log10 minimum covered CA sites per base",
            "matched-WT whole-cortex log2(CPM + 1)",
        ],
        "warning": (
            "P18 chromatin tracks are single source-study profiles and the methylomes have "
            "n=2 animals per genotype. Chromosome-block intervals measure genomic robustness, "
            "not animal-level uncertainty; associations are descriptive, not causal."
        ),
    }
    (output_dir / "mch_chromatin_context_models.json").write_text(
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
        default=Path("results/mch_gene_body_comparison/cortex_WT_expression_covariates.tsv"),
    )
    parser.add_argument(
        "--chromatin", type=Path,
        default=Path("results/mch_gene_body_comparison/cortex_chromatin_gene_signals.tsv.gz"),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/mch_gene_body_comparison")
    )
    parser.add_argument("--bootstraps", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=271828)
    args = parser.parse_args()
    print(json.dumps(build(args.methylation, args.expression, args.chromatin, args.output_dir, args.bootstraps, args.seed), indent=2))


if __name__ == "__main__":
    main()
