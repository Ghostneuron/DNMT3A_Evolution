#!/usr/bin/env python3
"""Test whether gene length predicts mCA effects beyond WT mCA and callability."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np


EXPRESSION_COLUMN = {
    "Dnmt3a1_KO-WT": "Dnmt3a1_KO_WT_mean_log2_CPM",
    "Dnmt3a2_KO-WT": "Dnmt3a2_KO_WT_mean_log2_CPM",
    "Dnmt3a1_delta_N-WT": "Dnmt3a1_delta_N_WT_mean_log2_CPM",
}


def zscore(values: np.ndarray) -> np.ndarray:
    scale = values.std(ddof=0)
    if scale == 0:
        raise ValueError("cannot standardize a constant variable")
    return (values - values.mean()) / scale


def ols(y: np.ndarray, predictors: list[np.ndarray]) -> tuple[np.ndarray, float]:
    design = np.column_stack([np.ones(len(y)), *predictors])
    beta, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ beta
    total = ((y - y.mean()) ** 2).sum()
    r_squared = 1.0 - float(residual @ residual) / float(total)
    return beta, r_squared


def percentile_interval(values: list[float]) -> tuple[float, float]:
    return tuple(float(x) for x in np.percentile(values, [2.5, 97.5]))


def analyze_contrast(
    rows: list[dict[str, str]], bootstraps: int, seed: int
) -> dict[str, object]:
    y = zscore(np.array([float(row["absolute_difference"]) for row in rows]))
    length = zscore(np.log10(np.array([float(row["length_bp"]) for row in rows])))
    lengths_bp = np.array([float(row["length_bp"]) for row in rows])
    minimum_sites = np.array(
        [float(row["minimum_sites_across_samples"]) for row in rows]
    )
    sites = zscore(np.log10(minimum_sites / lengths_bp))
    wt = zscore(np.array([float(row["wt_mean_corrected_mCA"]) for row in rows]))
    expression = zscore(np.array([float(row["matched_WT_log2_CPM"]) for row in rows]))
    chrom = np.array([row["chrom"] for row in rows])

    unadjusted, unadjusted_r2 = ols(y, [length])
    adjusted, adjusted_r2 = ols(y, [wt, length, sites, expression])
    reduced, reduced_r2 = ols(y, [wt, sites, expression])
    partial_r2 = max(0.0, (adjusted_r2 - reduced_r2) / (1.0 - reduced_r2))

    blocks = sorted(set(chrom))
    rng = np.random.default_rng(seed)
    boot_length = []
    for _ in range(bootstraps):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([np.flatnonzero(chrom == block) for block in selected])
        beta, _ = ols(
            y[indices],
            [wt[indices], length[indices], sites[indices], expression[indices]],
        )
        boot_length.append(float(beta[2]))
    ci_low, ci_high = percentile_interval(boot_length)

    return {
        "genes": len(rows),
        "chromosome_blocks": len(blocks),
        "unadjusted_standardized_length_beta": float(unadjusted[1]),
        "unadjusted_model_r_squared": unadjusted_r2,
        "adjusted_standardized_WT_mCA_beta": float(adjusted[1]),
        "adjusted_standardized_length_beta": float(adjusted[2]),
        "adjusted_standardized_CA_site_density_beta": float(adjusted[3]),
        "adjusted_standardized_expression_beta": float(adjusted[4]),
        "adjusted_model_r_squared": adjusted_r2,
        "length_partial_r_squared": partial_r2,
        "chromosome_block_bootstrap_length_beta_ci_low": ci_low,
        "chromosome_block_bootstrap_length_beta_ci_high": ci_high,
        "chromosome_block_bootstraps": bootstraps,
    }


def build(
    input_path: Path,
    expression_path: Path,
    output_dir: Path,
    bootstraps: int,
    seed: int,
) -> dict[str, object]:
    with expression_path.open(newline="") as handle:
        expression = {row["gene"]: row for row in csv.DictReader(handle, delimiter="\t")}
    grouped: dict[str, list[dict[str, str]]] = {}
    with gzip.open(input_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expression_row = expression.get(row["gene"])
            if expression_row is None:
                continue
            row["matched_WT_log2_CPM"] = expression_row[
                EXPRESSION_COLUMN[row["contrast"]]
            ]
            grouped.setdefault(row["contrast"], []).append(row)

    results = []
    for offset, contrast in enumerate(sorted(grouped)):
        result = analyze_contrast(grouped[contrast], bootstraps, seed + offset)
        results.append({"contrast": contrast, **result})

    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "gene_body_covariate_models.tsv"
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    report = {
        "complete": bool(results),
        "input": str(input_path),
        "expression_input": str(expression_path),
        "output": str(table_path),
        "outcome": "standardized genotype-minus-WT corrected gene-body mCA",
        "predictors": [
            "standardized WT corrected gene-body mCA",
            "standardized log10 gene length",
            "standardized log10 minimum covered CA sites per base across samples",
            "standardized matched-WT whole-cortex log2(CPM + 1)",
        ],
        "resampling": (
            "percentile interval from resampling autosomal chromosome blocks; "
            "this tests genomic robustness, not animal-level uncertainty"
        ),
        "warning": (
            "Observational secondary model with n=2 animals per genotype; coefficients "
            "must not be interpreted as causal or as independent biological replication."
        ),
        "results": results,
    }
    (output_dir / "gene_body_covariate_models.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results/mch_gene_body_comparison/gene_body_mCA_contrasts.tsv.gz"),
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("results/mch_gene_body_comparison/cortex_WT_expression_covariates.tsv"),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/mch_gene_body_comparison")
    )
    parser.add_argument("--bootstraps", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=314159)
    args = parser.parse_args()
    print(
        json.dumps(
            build(args.input, args.expression, args.output_dir, args.bootstraps, args.seed),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
