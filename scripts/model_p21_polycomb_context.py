#!/usr/bin/env python3
"""Test matched P21 Polycomb context, DNMT3A occupancy, and neuronal mCA effects."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np

try:
    from scripts.model_mch_gene_body_covariates import EXPRESSION_COLUMN, ols, zscore
except ModuleNotFoundError:
    from model_mch_gene_body_covariates import EXPRESSION_COLUMN, ols, zscore


MARK_FEATURES = (
    "h2ak119ub_gene_body_enrichment", "h2ak119ub_promoter_enrichment",
    "h3k27me3_gene_body_enrichment", "h3k27me3_promoter_enrichment",
    "h3k4me3_gene_body_enrichment", "h3k4me3_promoter_enrichment",
)


def analyze(
    outcome: np.ndarray,
    base: list[np.ndarray],
    focal: list[np.ndarray],
    focal_names: list[str],
    chromosomes: np.ndarray,
    bootstraps: int,
    seed: int,
) -> dict[str, object]:
    y = zscore(outcome)
    base_z = [zscore(value) for value in base]
    focal_z = [zscore(value) for value in focal]
    _, base_r2 = ols(y, base_z)
    extended, extended_r2 = ols(y, [*base_z, *focal_z])
    blocks = sorted(set(chromosomes))
    design = np.column_stack([np.ones(len(y)), *base_z, *focal_z])
    stats = {}
    for block in blocks:
        index = np.flatnonzero(chromosomes == block)
        x, block_y = design[index], y[index]
        stats[block] = (x.T @ x, x.T @ block_y, float(block_y @ block_y), float(block_y.sum()), len(block_y))
    rng = np.random.default_rng(seed)
    coefficient_samples = [[] for _ in focal]
    partial_samples = []
    base_columns = 1 + len(base)
    for _ in range(bootstraps):
        selected = rng.choice(blocks, len(blocks), replace=True)
        xtx = sum(stats[str(block)][0] for block in selected)
        xty = sum(stats[str(block)][1] for block in selected)
        yty = sum(stats[str(block)][2] for block in selected)
        ysum = sum(stats[str(block)][3] for block in selected)
        count = sum(stats[str(block)][4] for block in selected)
        beta = np.linalg.lstsq(xtx, xty, rcond=None)[0]
        total = max(1e-15, yty - ysum * ysum / count)
        full_r2 = 1.0 - max(0.0, yty - float(beta @ xty)) / total
        reduced_beta = np.linalg.lstsq(xtx[:base_columns, :base_columns], xty[:base_columns], rcond=None)[0]
        reduced_r2 = 1.0 - max(0.0, yty - float(reduced_beta @ xty[:base_columns])) / total
        for position, values in enumerate(coefficient_samples):
            values.append(float(beta[base_columns + position]))
        partial_samples.append(max(0.0, (full_r2 - reduced_r2) / max(1e-15, 1.0 - reduced_r2)))
    result: dict[str, object] = {
        "genes": len(y), "chromosome_blocks": len(blocks),
        "base_model_r_squared": base_r2, "extended_model_r_squared": extended_r2,
        "joint_partial_r_squared": max(0.0, (extended_r2 - base_r2) / max(1e-15, 1.0 - base_r2)),
        "joint_partial_r_squared_block_interval_low": float(np.percentile(partial_samples, 2.5)),
        "joint_partial_r_squared_block_interval_high": float(np.percentile(partial_samples, 97.5)),
        "chromosome_block_bootstraps": bootstraps,
    }
    for position, name in enumerate(focal_names):
        result[f"{name}_standardized_beta"] = float(extended[base_columns + position])
        result[f"{name}_beta_block_interval_low"] = float(np.percentile(coefficient_samples[position], 2.5))
        result[f"{name}_beta_block_interval_high"] = float(np.percentile(coefficient_samples[position], 97.5))
    return result


def add_features(signal: dict[str, str]) -> dict[str, float]:
    features = {}
    for mark in ("h2ak119ub", "h3k27me3", "h3k4me3"):
        for region in ("gene_body", "promoter"):
            features[f"{mark}_{region}_enrichment"] = float(
                np.log1p(float(signal[f"{mark}_{region}_mean"]))
                - np.log1p(float(signal[f"input_{region}_mean"]))
            )
    for region in ("gene_body", "promoter"):
        features[f"dnmt3a_occupancy_change_{region}"] = float(
            np.log1p(float(signal[f"delta_n_dnmt3a_{region}_mean"]))
            - np.log1p(float(signal[f"wt_dnmt3a_{region}_mean"]))
        )
        features[f"wt_dnmt3a_{region}"] = float(np.log1p(float(signal[f"wt_dnmt3a_{region}_mean"])))
    return features


def build(methylation_path: Path, expression_path: Path, signal_path: Path, output_dir: Path, bootstraps: int, seed: int) -> dict[str, object]:
    with expression_path.open(newline="") as handle:
        expression = {row["gene"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with gzip.open(signal_path, "rt", newline="") as handle:
        signals = {(row["chrom"], row["start"], row["end"], row["gene"]): add_features(row) for row in csv.DictReader(handle, delimiter="\t")}
    grouped: dict[str, list[dict[str, object]]] = {}
    with gzip.open(methylation_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expr = expression.get(row["gene"])
            features = signals.get((row["chrom"], row["start"], row["end"], row["gene"]))
            if expr is None or features is None:
                continue
            row["matched_WT_log2_CPM"] = expr[EXPRESSION_COLUMN[row["contrast"]]]
            grouped.setdefault(row["contrast"], []).append({**row, **features})

    single_results = []
    joint_results = []
    iteration = 0
    for contrast in sorted(grouped):
        for sensitivity, rows in (
            ("all_genes", grouped[contrast]),
            ("WT_mCA_ge_0.005", [row for row in grouped[contrast] if float(row["wt_mean_corrected_mCA"]) >= 0.005]),
        ):
            lengths = np.array([float(row["length_bp"]) for row in rows])
            base = [
                np.array([float(row["wt_mean_corrected_mCA"]) for row in rows]),
                np.log10(lengths),
                np.log10(np.array([float(row["minimum_sites_across_samples"]) for row in rows]) / lengths),
                np.array([float(row["matched_WT_log2_CPM"]) for row in rows]),
            ]
            outcome = np.array([float(row["absolute_difference"]) for row in rows])
            chromosomes = np.array([str(row["chrom"]) for row in rows])
            for feature in (*MARK_FEATURES, "dnmt3a_occupancy_change_gene_body", "dnmt3a_occupancy_change_promoter"):
                focal = np.array([float(row[feature]) for row in rows])
                result = analyze(outcome, base, [focal], [feature], chromosomes, bootstraps, seed + iteration)
                single_results.append({"contrast": contrast, "sensitivity": sensitivity, "feature": feature, **result})
                iteration += 1
            for region in ("gene_body", "promoter"):
                names = [f"h2ak119ub_{region}_enrichment", f"h3k27me3_{region}_enrichment"]
                result = analyze(outcome, base, [np.array([float(row[name]) for row in rows]) for name in names], names, chromosomes, bootstraps, seed + iteration)
                joint_results.append({"contrast": contrast, "sensitivity": sensitivity, "region": region, "model": "h2ak119ub_plus_h3k27me3", **result})
                iteration += 1
                names = [f"h2ak119ub_{region}_enrichment", f"h3k27me3_{region}_enrichment", f"h3k4me3_{region}_enrichment"]
                result = analyze(outcome, base, [np.array([float(row[name]) for row in rows]) for name in names], names, chromosomes, bootstraps, seed + iteration)
                joint_results.append({"contrast": contrast, "sensitivity": sensitivity, "region": region, "model": "three_marks", **result})
                iteration += 1

    occupancy_rows = grouped["Dnmt3a1_delta_N-WT"]
    occupancy_results = []
    lengths = np.array([float(row["length_bp"]) for row in occupancy_rows])
    chromosomes = np.array([str(row["chrom"]) for row in occupancy_rows])
    for region in ("gene_body", "promoter"):
        outcome = np.array([float(row[f"dnmt3a_occupancy_change_{region}"]) for row in occupancy_rows])
        base = [
            np.array([float(row[f"wt_dnmt3a_{region}"]) for row in occupancy_rows]),
            np.log10(lengths),
            np.array([float(row["matched_WT_log2_CPM"]) for row in occupancy_rows]),
        ]
        for mark in ("h2ak119ub", "h3k27me3", "h3k4me3"):
            name = f"{mark}_{region}_enrichment"
            result = analyze(outcome, base, [np.array([float(row[name]) for row in occupancy_rows])], [name], chromosomes, bootstraps, seed + iteration)
            occupancy_results.append({"region": region, "model": mark, **result})
            iteration += 1
        names = [f"h2ak119ub_{region}_enrichment", f"h3k27me3_{region}_enrichment"]
        result = analyze(outcome, base, [np.array([float(row[name]) for row in occupancy_rows]) for name in names], names, chromosomes, bootstraps, seed + iteration)
        occupancy_results.append({"region": region, "model": "h2ak119ub_plus_h3k27me3", **result})
        iteration += 1
        names = [f"h2ak119ub_{region}_enrichment", f"h3k27me3_{region}_enrichment", f"h3k4me3_{region}_enrichment"]
        result = analyze(outcome, base, [np.array([float(row[name]) for row in occupancy_rows]) for name in names], names, chromosomes, bootstraps, seed + iteration)
        occupancy_results.append({"region": region, "model": "three_marks", **result})
        iteration += 1

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for filename, rows in (
        ("p21_polycomb_mca_single_models.tsv", single_results),
        ("p21_polycomb_mca_joint_models.tsv", joint_results),
        ("p21_dnmt3a_occupancy_redistribution_models.tsv", occupancy_results),
    ):
        path = output_dir / filename
        fields = []
        for row in rows:
            for field in row:
                if field not in fields:
                    fields.append(field)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader(); writer.writerows(rows)
        outputs.append(str(path))
    correlation_path = output_dir / "p21_polycomb_feature_correlations.tsv"
    correlation_features = list(MARK_FEATURES)
    with correlation_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["feature_1", "feature_2", "pearson_r"], delimiter="\t")
        writer.writeheader()
        for first_index, first in enumerate(correlation_features):
            first_values = np.array([float(row[first]) for row in occupancy_rows])
            for second in correlation_features[first_index:]:
                second_values = np.array([float(row[second]) for row in occupancy_rows])
                writer.writerow({"feature_1": first, "feature_2": second, "pearson_r": float(np.corrcoef(first_values, second_values)[0, 1])})
    outputs.append(str(correlation_path))
    report = {
        "complete": True, "outputs": outputs,
        "mca_single_models": len(single_results), "mca_joint_models": len(joint_results),
        "occupancy_redistribution_models": len(occupancy_results),
        "interpretation_boundary": (
            "P21 neuronal histone marks are age/cell matched to the methylomes, but each bigWig is a single source-study profile. "
            "DNMT3A occupancy is from whole cortex and has one track per genotype. Chromosome-block intervals measure genomic, not animal-level, robustness."
        ),
    }
    (output_dir / "p21_polycomb_models.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--methylation", type=Path, default=Path("results/mch_gene_body_comparison/gene_body_mCA_contrasts.tsv.gz"))
    parser.add_argument("--expression", type=Path, default=Path("results/mch_gene_body_comparison/cortex_WT_expression_covariates.tsv"))
    parser.add_argument("--signals", type=Path, default=Path("results/mch_gene_body_comparison/p21_polycomb_gene_signals.tsv.gz"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/mch_gene_body_comparison"))
    parser.add_argument("--bootstraps", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=141421)
    args = parser.parse_args()
    print(json.dumps(build(args.methylation, args.expression, args.signals, args.output_dir, args.bootstraps, args.seed), indent=2))


if __name__ == "__main__":
    main()
