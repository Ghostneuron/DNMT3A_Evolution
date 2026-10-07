#!/usr/bin/env python3
"""Test whether the dunnart P12/P20 contrast depends on one animal."""

from __future__ import annotations

import csv
import json
import xml.etree.ElementTree as ET
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/06_isoform_evolution"
JUNCTIONS = RESULTS / "dunnart_neocortex_stage_junctions.tsv"
SALMON = RESULTS / "dunnart_salmon_replicate_metrics.tsv"
SRA_XML = ROOT / "data/raw/marsupial_sra/dunnart_sra_records.xml"
RUNS = {f"SRR1303604{index}" for index in range(3, 9)}
METRICS = (
    "junction_internal_rate_fold",
    "junction_internal_common_ratio_fold",
    "junction_internal_full_ratio_fold",
    "salmon_internal_rate_fold",
    "salmon_internal_fraction_fold",
    "salmon_internal_full_ratio_fold",
)


def number(row: dict, key: str) -> float:
    return float(row[key])


def pooled_metrics(rows: list[dict]) -> dict[str, float]:
    pairs = sum(number(row, "paired_read_count") for row in rows)
    internal_pairs = sum(
        number(row, "internal_067_first_junction_pairs") for row in rows
    )
    common_pairs = sum(
        number(row, "common_downstream_junction_pairs") for row in rows
    )
    full_pairs = sum(
        number(row, "full_length_first_junction_pairs") for row in rows
    )
    input_fragments = sum(number(row, "input_fragments") for row in rows)
    retained = sum(number(row, "retained_dnmt3a_fragments") for row in rows)
    internal_fragments = sum(
        number(row, "internal_066_067_estimated_fragments") for row in rows
    )
    full_fragments = sum(
        number(row, "internal_066_067_estimated_fragments")
        / number(row, "internal_066_067_to_full_length_count_ratio")
        for row in rows
    )
    return {
        "junction_internal_rate": internal_pairs / pairs,
        "junction_internal_common_ratio": internal_pairs / common_pairs,
        "junction_internal_full_ratio": internal_pairs / full_pairs,
        "salmon_internal_rate": internal_fragments / input_fragments,
        "salmon_internal_fraction": internal_fragments / retained,
        "salmon_internal_full_ratio": internal_fragments / full_fragments,
    }


def contrast(rows: list[dict]) -> dict[str, float]:
    stage_metrics = {
        stage: pooled_metrics([row for row in rows if row["stage"] == stage])
        for stage in ("P12", "P20")
    }
    return {
        f"{metric}_fold": stage_metrics["P12"][metric] / stage_metrics["P20"][metric]
        for metric in stage_metrics["P12"]
    }


def metadata_rows() -> list[dict[str, str]]:
    root = ET.parse(SRA_XML).getroot()
    rows = []
    for package in root.findall("EXPERIMENT_PACKAGE"):
        experiment = package.find("EXPERIMENT")
        sample = package.find("SAMPLE")
        if experiment is None or sample is None:
            continue
        attributes = {}
        for item in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
            key = (item.findtext("TAG") or "").strip().lower()
            value = (item.findtext("VALUE") or "").strip()
            attributes[key] = value
        descriptor = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR")
        platform = experiment.find("PLATFORM")
        platform_name = next(iter(platform)).tag if platform is not None else ""
        for run in package.findall("./RUN_SET/RUN"):
            accession = run.get("accession", "")
            if accession not in RUNS:
                continue
            rows.append({
                "run_accession": accession,
                "sample_title": sample.findtext("TITLE") or "",
                "sample_accession": sample.get("accession", ""),
                "tissue": attributes.get("tissue", ""),
                "developmental_stage": attributes.get(
                    "developmental stage", ""
                ),
                "sex": attributes.get("sex", ""),
                "genotype": attributes.get("genotype", ""),
                "library_strategy": (
                    descriptor.findtext("LIBRARY_STRATEGY")
                    if descriptor is not None else ""
                ),
                "library_selection": (
                    descriptor.findtext("LIBRARY_SELECTION")
                    if descriptor is not None else ""
                ),
                "library_layout": (
                    next(iter(descriptor.find("LIBRARY_LAYOUT"))).tag
                    if descriptor is not None
                    and descriptor.find("LIBRARY_LAYOUT") is not None else ""
                ),
                "platform": platform_name,
                "instrument": (
                    platform.findtext(f"{platform_name}/INSTRUMENT_MODEL")
                    if platform is not None else ""
                ),
            })
    return sorted(rows, key=lambda row: row["run_accession"])


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    with JUNCTIONS.open() as handle:
        junction_rows = {
            row["run_accession"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    with SALMON.open() as handle:
        salmon_rows = {
            row["run_accession"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    if set(junction_rows) != RUNS or set(salmon_rows) != RUNS:
        raise SystemExit("Robustness analysis requires the complete six-run set")
    rows = []
    for run in sorted(RUNS):
        rows.append({**junction_rows[run], **salmon_rows[run]})

    observed = contrast(rows)
    single_deletions = []
    for omitted in sorted(RUNS):
        values = contrast([row for row in rows if row["run_accession"] != omitted])
        single_deletions.append({
            "omitted_run": omitted,
            "omitted_stage": next(
                row["stage"] for row in rows if row["run_accession"] == omitted
            ),
            **{key: round(value, 6) for key, value in values.items()},
        })

    p12_runs = sorted(row["run_accession"] for row in rows if row["stage"] == "P12")
    p20_runs = sorted(row["run_accession"] for row in rows if row["stage"] == "P20")
    balanced_deletions = []
    for omitted_p12, omitted_p20 in product(p12_runs, p20_runs):
        values = contrast([
            row for row in rows
            if row["run_accession"] not in {omitted_p12, omitted_p20}
        ])
        balanced_deletions.append({
            "omitted_P12_run": omitted_p12,
            "omitted_P20_run": omitted_p20,
            **{key: round(value, 6) for key, value in values.items()},
        })

    metadata = metadata_rows()
    technical_fields = (
        "tissue", "genotype", "library_strategy", "library_selection",
        "library_layout", "platform", "instrument",
    )
    summary = {
        "observed_three_by_three_fold_contrasts": {
            key: round(value, 6) for key, value in observed.items()
        },
        "single_animal_deletion": {
            "scenario_count": len(single_deletions),
            "fold_ranges": {
                metric: [
                    min(row[metric] for row in single_deletions),
                    max(row[metric] for row in single_deletions),
                ]
                for metric in METRICS
            },
            "all_contrasts_above_one": {
                metric: all(row[metric] > 1 for row in single_deletions)
                for metric in METRICS
            },
        },
        "balanced_one_per_stage_deletion": {
            "scenario_count": len(balanced_deletions),
            "fold_ranges": {
                metric: [
                    min(row[metric] for row in balanced_deletions),
                    max(row[metric] for row in balanced_deletions),
                ]
                for metric in METRICS
            },
            "all_contrasts_above_one": {
                metric: all(row[metric] > 1 for row in balanced_deletions)
                for metric in METRICS
            },
        },
        "metadata_audit": {
            "distinct_sample_accessions": len({
                row["sample_accession"] for row in metadata
            }),
            "technical_fields_uniform_across_all_six": {
                field: len({row[field] for row in metadata}) == 1
                for field in technical_fields
            },
            "missing_sex_records": sum(not row["sex"] for row in metadata),
            "sex_confounding_excluded": all(row["sex"] for row in metadata),
            "explicit_library_batch_field_available": False,
        },
        "interpretation": (
            "Direction is robust to every single-animal deletion and every "
            "balanced one-per-stage deletion only for metrics whose "
            "all_contrasts_above_one flag is true. Missing sex and library-batch "
            "metadata prevent exclusion of those stage-correlated confounders."
        ),
    }
    write_tsv(
        RESULTS / "dunnart_neocortex_single_animal_deletion.tsv",
        single_deletions,
    )
    write_tsv(
        RESULTS / "dunnart_neocortex_balanced_deletion.tsv",
        balanced_deletions,
    )
    write_tsv(RESULTS / "dunnart_neocortex_sample_metadata.tsv", metadata)
    (
        RESULTS / "dunnart_neocortex_stage_robustness_summary.json"
    ).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
