#!/usr/bin/env python3
"""Localize replicated DNMT3A genotype effects in gene-body mCA summaries.

This is a descriptive secondary analysis. It deliberately reports no
cytosite- or gene-level p-values because genomic observations do not replace
the two biological replicates in each genotype.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
from collections import defaultdict
from pathlib import Path

try:
    from scripts.run_mch_bismark import DEFAULT_EXTERNAL_ROOT, read_samples
except ModuleNotFoundError:
    from run_mch_bismark import DEFAULT_EXTERNAL_ROOT, read_samples


GROUPS = ("WT", "Dnmt3a1_KO", "Dnmt3a2_KO", "Dnmt3a1_delta_N")
AUTOSOMES = {f"chr{number}" for number in range(1, 20)}


def corrected(fraction: float, background: float) -> float:
    if background >= 1:
        raise ValueError("background must be below one")
    return max(0.0, (fraction - background) / (1.0 - background))


def median_or_blank(values: list[float]) -> float | str:
    return statistics.median(values) if values else ""


def lambda_ca_background(path: Path) -> float:
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["scope"] == "phage_lambda" and row["context"] == "CA":
                return float(row["fraction"])
    raise ValueError(f"no lambda CA row in {path}")


def read_gene_ca(path: Path) -> dict[tuple[str, int, int, str, str], tuple[int, float]]:
    rows = {}
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["context"] != "CA" or row["chrom"] not in AUTOSOMES:
                continue
            key = (
                row["chrom"],
                int(row["start"]),
                int(row["end"]),
                row["gene"],
                row["strand"],
            )
            rows[key] = (int(row["sites"]), float(row["fraction"]))
    return rows


def build(
    external_root: Path,
    output_dir: Path,
    min_sites: int = 100,
) -> dict[str, object]:
    samples = read_samples()
    values: dict[tuple[str, int, int, str, str], dict[str, tuple[int, float]]] = (
        defaultdict(dict)
    )
    backgrounds = {}
    for run in samples:
        source = external_root / "results/bismark" / run / "context_summary"
        backgrounds[run] = lambda_ca_background(source / "context_summary.tsv")
        for key, (sites, fraction) in read_gene_ca(
            source / "gene_body_context_summary.tsv.gz"
        ).items():
            values[key][run] = (sites, corrected(fraction, backgrounds[run]))

    expected = set(samples)
    rows = []
    for key, run_values in values.items():
        if set(run_values) != expected:
            continue
        if min(sites for sites, _ in run_values.values()) < min_sites:
            continue
        chrom, start, end, gene, strand = key
        group_values = {
            group: [
                run_values[run][1]
                for run, sample in samples.items()
                if sample["genotype"] == group
            ]
            for group in GROUPS
        }
        wt = statistics.mean(group_values["WT"])
        base = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "gene": gene,
            "strand": strand,
            "length_bp": end - start,
            "long_gene_ge_100kb": end - start >= 100_000,
            "minimum_sites_across_samples": min(
                sites for sites, _ in run_values.values()
            ),
            "wt_mean_corrected_mCA": wt,
        }
        for group in GROUPS[1:]:
            treatment = statistics.mean(group_values[group])
            rows.append(
                {
                    **base,
                    "contrast": f"{group}-WT",
                    "treatment_mean_corrected_mCA": treatment,
                    "absolute_difference": treatment - wt,
                    "relative_to_WT": treatment / wt if wt > 0 else "",
                    "wt_replicate_1": group_values["WT"][0],
                    "wt_replicate_2": group_values["WT"][1],
                    "treatment_replicate_1": group_values[group][0],
                    "treatment_replicate_2": group_values[group][1],
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "gene_body_mCA_contrasts.tsv.gz"
    fields = list(rows[0]) if rows else []
    with gzip.open(table_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary_rows = []
    for contrast in (f"{group}-WT" for group in GROUPS[1:]):
        selected = [row for row in rows if row["contrast"] == contrast]
        effects = [float(row["absolute_difference"]) for row in selected]
        ratios = [
            float(row["relative_to_WT"])
            for row in selected
            if row["relative_to_WT"] != ""
        ]
        long_selected = [row for row in selected if row["long_gene_ge_100kb"]]
        short_selected = [row for row in selected if not row["long_gene_ge_100kb"]]
        long_ratios = [
            float(row["relative_to_WT"])
            for row in long_selected
            if row["relative_to_WT"] != ""
        ]
        short_ratios = [
            float(row["relative_to_WT"])
            for row in short_selected
            if row["relative_to_WT"] != ""
        ]
        long_effects = [float(row["absolute_difference"]) for row in long_selected]
        short_effects = [float(row["absolute_difference"]) for row in short_selected]
        summary_rows.append(
            {
                "contrast": contrast,
                "genes": len(selected),
                "median_absolute_difference": statistics.median(effects),
                "median_relative_to_WT": statistics.median(ratios),
                "genes_ge_100kb": sum(
                    bool(row["long_gene_ge_100kb"]) for row in selected
                ),
                "median_relative_to_WT_genes_ge_100kb": median_or_blank(long_ratios),
                "median_relative_to_WT_genes_lt_100kb": median_or_blank(short_ratios),
                "median_absolute_difference_genes_ge_100kb": median_or_blank(
                    long_effects
                ),
                "median_absolute_difference_genes_lt_100kb": median_or_blank(
                    short_effects
                ),
            }
        )
    summary_path = output_dir / "gene_body_mCA_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(summary_rows[0]) if summary_rows else [], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    report = {
        "complete": len(samples) == 8,
        "samples": len(samples),
        "minimum_CA_sites_per_gene_per_sample": min_sites,
        "complete_case_genes": len(rows) // 3,
        "contrast_table": str(table_path),
        "summary_table": str(summary_path),
        "warning": (
            "Descriptive secondary analysis; genes and cytosites are not "
            "independent biological replicates, so no gene-level p-values are reported."
        ),
    }
    (output_dir / "gene_body_mCA_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/mch_gene_body_comparison")
    )
    parser.add_argument("--min-sites", type=int, default=100)
    args = parser.parse_args()
    print(json.dumps(build(args.external_root, args.output_dir, args.min_sites), indent=2))


if __name__ == "__main__":
    main()
