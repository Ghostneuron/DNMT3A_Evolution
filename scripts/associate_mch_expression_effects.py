#!/usr/bin/env python3
"""Associate gene-body mCA effects with matched cortex expression effects."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

try:
    from scripts.model_mch_gene_body_covariates import ols, percentile_interval, zscore
except ModuleNotFoundError:
    from model_mch_gene_body_covariates import ols, percentile_interval, zscore


def analyze(rows: list[dict[str, str]], bootstraps: int, seed: int) -> dict[str, object]:
    expression_effect = zscore(
        np.array([float(row["expression_difference"]) for row in rows])
    )
    mca_effect = zscore(np.array([float(row["absolute_difference"]) for row in rows]))
    wt_expression = zscore(np.array([float(row["wt_mean_log2_CPM"]) for row in rows]))
    wt_mca = zscore(np.array([float(row["wt_mean_corrected_mCA"]) for row in rows]))
    length_bp = np.array([float(row["length_bp"]) for row in rows])
    length = zscore(np.log10(length_bp))
    density = zscore(
        np.log10(
            np.array([float(row["minimum_sites_across_samples"]) for row in rows])
            / length_bp
        )
    )
    chrom = np.array([row["chrom"] for row in rows])

    unadjusted_rho = float(spearmanr(mca_effect, expression_effect).statistic)
    full, full_r2 = ols(
        expression_effect, [mca_effect, wt_expression, wt_mca, length, density]
    )
    reduced, reduced_r2 = ols(expression_effect, [wt_expression, wt_mca, length, density])
    partial_r2 = max(0.0, (full_r2 - reduced_r2) / (1.0 - reduced_r2))

    blocks = sorted(set(chrom))
    rng = np.random.default_rng(seed)
    coefficients = []
    for _ in range(bootstraps):
        selected = rng.choice(blocks, len(blocks), replace=True)
        indices = np.concatenate([np.flatnonzero(chrom == block) for block in selected])
        beta, _ = ols(
            expression_effect[indices],
            [
                mca_effect[indices],
                wt_expression[indices],
                wt_mca[indices],
                length[indices],
                density[indices],
            ],
        )
        coefficients.append(float(beta[1]))
    low, high = percentile_interval(coefficients)
    return {
        "genes": len(rows),
        "chromosome_blocks": len(blocks),
        "unadjusted_spearman_rho": unadjusted_rho,
        "adjusted_standardized_mCA_effect_beta": float(full[1]),
        "adjusted_model_r_squared": full_r2,
        "mCA_effect_partial_r_squared": partial_r2,
        "chromosome_block_beta_ci_low": low,
        "chromosome_block_beta_ci_high": high,
        "chromosome_block_bootstraps": bootstraps,
    }


def build(
    methylation_path: Path,
    expression_path: Path,
    output_dir: Path,
    bootstraps: int,
    seed: int,
) -> dict[str, object]:
    expression = {}
    with gzip.open(expression_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expression[(row["contrast"], row["gene"])] = row
    grouped = {}
    with gzip.open(methylation_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            effect = expression.get((row["contrast"], row["gene"]))
            if effect is None:
                continue
            merged = {**row, **effect}
            grouped.setdefault(row["contrast"], []).append(merged)

    results = [
        {"contrast": contrast, **analyze(grouped[contrast], bootstraps, seed + offset)}
        for offset, contrast in enumerate(sorted(grouped))
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "mCA_expression_effect_associations.tsv"
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    report = {
        "complete": bool(results),
        "outcome": "standardized matched-cortex expression difference",
        "focal_predictor": "standardized genotype-minus-WT corrected gene-body mCA",
        "covariates": ["WT expression", "WT mCA", "log gene length", "covered CA-site density"],
        "output": str(table_path),
        "warning": (
            "Cross-assay descriptive association. Whole-cortex RNA and NeuN-positive "
            "methylomes are not cell- or animal-matched; n=2 methylomes per genotype."
        ),
        "results": results,
    }
    (output_dir / "mCA_expression_effect_associations.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--methylation",
        type=Path,
        default=Path("results/mch_gene_body_comparison/gene_body_mCA_contrasts.tsv.gz"),
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("results/mch_gene_body_comparison/cortex_expression_contrasts.tsv.gz"),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/mch_gene_body_comparison")
    )
    parser.add_argument("--bootstraps", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=271828)
    args = parser.parse_args()
    print(
        json.dumps(
            build(args.methylation, args.expression, args.output_dir, args.bootstraps, args.seed),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
