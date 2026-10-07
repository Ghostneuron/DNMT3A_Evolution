#!/usr/bin/env python3
"""Compare context-resolved methylation summaries across DNMT3A genotypes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

try:
    from scripts.run_mch_bismark import DEFAULT_EXTERNAL_ROOT, read_samples
except ModuleNotFoundError:  # Direct execution from the project root.
    from run_mch_bismark import DEFAULT_EXTERNAL_ROOT, read_samples


def read_sample_summary(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(newline="") as handle:
        return {
            (row["scope"], row["context"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def mean(values: list[float]) -> float:
    return statistics.mean(values) if values else float("nan")


def sd(values: list[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else float("nan")


def compare(
    external_root: Path,
    summaries: dict[str, Path],
    output_dir: Path,
) -> dict[str, object]:
    samples = read_samples()
    sample_rows = []
    for run, path in summaries.items():
        table = read_sample_summary(path)
        sample = samples[run]
        for (scope, context), row in table.items():
            if scope not in {"autosome", "genome", "phage_lambda", "plasmid_puc19c"}:
                continue
            fraction = float(row["fraction"])
            lambda_row = table.get(("phage_lambda", context)) or table.get(
                ("phage_lambda", "all_C")
            )
            conversion_error = (
                float(lambda_row["fraction"]) if lambda_row is not None else math.nan
            )
            is_non_cpg = context in {"CA", "CC", "CT", "CH", "CHG", "CHH"}
            corrected = (
                max(0.0, (fraction - conversion_error) / (1.0 - conversion_error))
                if is_non_cpg
                and scope in {"autosome", "genome"}
                and not math.isnan(conversion_error)
                and conversion_error < 1.0
                else fraction
            )
            sample_rows.append(
                {
                    "run_accession": run,
                    "sample_accession": sample["sample_accession"],
                    "genotype": sample["genotype"],
                    "replicate": sample["replicate"],
                    "scope": scope,
                    "context": context,
                    "methylated": int(row["methylated"]),
                    "unmethylated": int(row["unmethylated"]),
                    "sites": int(row["sites"]),
                    "fraction": fraction,
                    "conversion_error": conversion_error,
                    "analysis_fraction": corrected,
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    sample_path = output_dir / "sample_context_summary.tsv"
    sample_fields = list(sample_rows[0]) if sample_rows else []
    with sample_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_rows)

    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in sample_rows:
        grouped[(row["genotype"], row["scope"], row["context"])].append(row)

    group_rows = []
    for (genotype, scope, context), rows in sorted(grouped.items()):
        fractions = [float(row["analysis_fraction"]) for row in rows]
        group_rows.append(
            {
                "genotype": genotype,
                "scope": scope,
                "context": context,
                "n": len(rows),
                "mean_fraction": f"{mean(fractions):.10g}",
                "sd_fraction": f"{sd(fractions):.10g}",
                "min_fraction": f"{min(fractions):.10g}",
                "max_fraction": f"{max(fractions):.10g}",
            }
        )
    group_path = output_dir / "genotype_context_summary.tsv"
    with group_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(group_rows[0]) if group_rows else [], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(group_rows)

    lookup = {
        (row["genotype"], row["scope"], row["context"]): row
        for row in group_rows
    }
    contrast_rows = []
    for key, treatment in sorted(lookup.items()):
        genotype, scope, context = key
        if genotype == "WT":
            continue
        control = lookup.get(("WT", scope, context))
        if control is None:
            continue
        effect = float(treatment["mean_fraction"]) - float(
            control["mean_fraction"]
        )
        contrast_rows.append(
            {
                "contrast": f"{genotype}-WT",
                "scope": scope,
                "context": context,
                "treatment_n": treatment["n"],
                "wt_n": control["n"],
                "treatment_mean": treatment["mean_fraction"],
                "wt_mean": control["mean_fraction"],
                "absolute_difference": f"{effect:.10g}",
                "relative_difference": f"{(
                    effect / float(control['mean_fraction'])
                    if float(control['mean_fraction']) != 0
                    else math.nan
                ):.10g}",
            }
        )
    contrast_path = output_dir / "genotype_vs_wt_contrasts.tsv"
    with contrast_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(contrast_rows[0]) if contrast_rows else [],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(contrast_rows)

    report = {
        "external_root": str(external_root),
        "samples_found": len(summaries),
        "samples_expected": len(samples),
        "complete": len(summaries) == len(samples),
        "sample_table": str(sample_path),
        "group_table": str(group_path),
        "contrast_table": str(contrast_path),
    }
    (output_dir / "comparison_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def discover_summaries(external_root: Path) -> dict[str, Path]:
    summaries = {}
    for run in read_samples():
        candidates = list(
            (
                external_root
                / "results/bismark"
                / run
                / "context_summary"
            ).glob("context_summary.tsv")
        )
        if candidates:
            summaries[run] = candidates[0]
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output = args.output_dir or (
        args.external_root / "results/bismark/genotype_comparison"
    )
    print(
        json.dumps(
            compare(
                args.external_root,
                discover_summaries(args.external_root),
                output,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
