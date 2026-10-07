#!/usr/bin/env python3
"""Build final replicate-level QC and genotype effect summaries for EM-seq."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

try:
    from scripts.run_mch_bismark import DEFAULT_EXTERNAL_ROOT, read_samples
except ModuleNotFoundError:
    from run_mch_bismark import DEFAULT_EXTERNAL_ROOT, read_samples


GROUP_ORDER = ("WT", "Dnmt3a1_KO", "Dnmt3a2_KO", "Dnmt3a1_delta_N")


def match_number(text: str, pattern: str, cast: type = float) -> float | int:
    match = re.search(pattern, text, flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"pattern not found: {pattern}")
    return cast(match.group(1))


def read_context_rows(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    with path.open(newline="") as handle:
        return {
            (row["run_accession"], row["scope"], row["context"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def mean(values: list[float]) -> float:
    return statistics.mean(values)


def sd(values: list[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else math.nan


def build(external_root: Path, output_dir: Path) -> dict[str, object]:
    samples = read_samples()
    comparison = external_root / "results/bismark/genotype_comparison"
    contexts = read_context_rows(comparison / "sample_context_summary.tsv")
    if len({key[0] for key in contexts}) != len(samples):
        raise RuntimeError("final context comparison does not contain all samples")

    rows: list[dict[str, object]] = []
    for run, sample in samples.items():
        root = external_root / "results/bismark" / run
        alignment = (root / "alignment" / f"{run}_PE_report.txt").read_text()
        dedup = (
            root / "deduplicated" / f"{run}_pe.deduplication_report.txt"
        ).read_text()
        splitting = (
            root
            / "methylation"
            / f"{run}_pe.deduplicated_splitting_report.txt"
        ).read_text()
        ca = contexts[(run, "autosome", "CA")]
        ch = contexts[(run, "autosome", "CH")]
        rows.append(
            {
                "run_accession": run,
                "genotype": sample["genotype"],
                "replicate": int(sample["replicate"]),
                "input_pairs": match_number(
                    alignment, r"Sequence pairs analysed in total:\s*(\d+)", int
                ),
                "mapping_efficiency": match_number(
                    alignment, r"Mapping efficiency:\s*([0-9.]+)%"
                )
                / 100,
                "deduplicated_pairs": match_number(
                    dedup, r"deduplicated leftover sequences:\s*(\d+)", int
                ),
                "duplicate_fraction": match_number(
                    dedup, r"removed:\s*\d+\s*\(([0-9.]+)%\)"
                )
                / 100,
                "mCG": match_number(
                    splitting, r"C methylated in CpG context:\s*([0-9.]+)%"
                )
                / 100,
                "mCHG": match_number(
                    splitting, r"C methylated in CHG context:\s*([0-9.]+)%"
                )
                / 100,
                "mCHH": match_number(
                    splitting, r"C methylated in CHH context:\s*([0-9.]+)%"
                )
                / 100,
                "autosomal_raw_mCA": float(ca["fraction"]),
                "lambda_mCA_background": float(ca["conversion_error"]),
                "autosomal_corrected_mCA": float(ca["analysis_fraction"]),
                "autosomal_raw_mCH": float(ch["fraction"]),
                "lambda_mCH_background": float(ch["conversion_error"]),
                "autosomal_corrected_mCH": float(ch["analysis_fraction"]),
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    replicate_path = output_dir / "replicate_qc_and_methylation.tsv"
    with replicate_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["genotype"])].append(row)
    wt_mca = mean([float(row["autosomal_corrected_mCA"]) for row in grouped["WT"]])
    wt_mch = mean([float(row["autosomal_corrected_mCH"]) for row in grouped["WT"]])
    group_rows: list[dict[str, object]] = []
    for genotype in GROUP_ORDER:
        members = grouped[genotype]
        mca = [float(row["autosomal_corrected_mCA"]) for row in members]
        mch = [float(row["autosomal_corrected_mCH"]) for row in members]
        mca_mean = mean(mca)
        mch_mean = mean(mch)
        group_rows.append(
            {
                "genotype": genotype,
                "n": len(members),
                "corrected_mCA_mean": mca_mean,
                "corrected_mCA_sd": sd(mca),
                "corrected_mCA_difference_vs_WT": mca_mean - wt_mca,
                "corrected_mCA_relative_vs_WT": mca_mean / wt_mca,
                "corrected_mCH_mean": mch_mean,
                "corrected_mCH_sd": sd(mch),
                "corrected_mCH_difference_vs_WT": mch_mean - wt_mch,
                "corrected_mCH_relative_vs_WT": mch_mean / wt_mch,
            }
        )
    group_path = output_dir / "genotype_effect_summary.tsv"
    with group_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(group_rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(group_rows)

    report = {
        "complete": len(rows) == len(samples) == 8,
        "samples": len(rows),
        "replicates_per_genotype": 2,
        "primary_endpoint": "lambda-corrected autosomal mCA fraction",
        "replicate_table": str(replicate_path),
        "group_table": str(group_path),
        "interpretation": (
            "Effect sizes are descriptive because n=2 per genotype. A value clamped "
            "to zero means no signal above the matched lambda estimate, not proof of "
            "literal absence of biological methylation."
        ),
    }
    (output_dir / "final_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output = args.output_dir or (
        args.external_root / "results/bismark/final_summary"
    )
    print(json.dumps(build(args.external_root, output), indent=2))


if __name__ == "__main__":
    main()
