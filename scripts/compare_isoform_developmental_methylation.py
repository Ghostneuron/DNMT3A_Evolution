#!/usr/bin/env python3
"""Compare Dnmt3a1/2 knockout methylation effects across development.

Uses the processed MM285 beta-value matrix from GSE295720. Contrasts are
estimated within tissue and developmental stage. Probe-level values summarize
effect distributions; probes are not treated as biological replicates.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import statistics
from collections import defaultdict
from pathlib import Path

try:
    from scripts.audit_isoform_functional_datasets import parse_soft_samples
except ModuleNotFoundError:  # direct execution from the scripts directory
    from audit_isoform_functional_datasets import parse_soft_samples


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/isoform_functional_genomics/GSE295720"
MATRIX = RAW / "GSE295720_MM285_betavalue_all.txt.gz"
METADATA = RAW / "GSE295720_samples.soft.txt"
SOURCE_MANIFEST = RAW / "source_manifest.tsv"
RESULTS = ROOT / "results/07_isoform_function/GSE295720_developmental_methylation"
STAGES = ("E15.5", "PD21")
TISSUES = ("brain", "liver")
KNOCKOUTS = ("Dnmt3a1-/-", "Dnmt3a2-/-")
EFFECT_THRESHOLD = 0.10


def validate_sources() -> None:
    with SOURCE_MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        path = RAW / row["local_file"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != row["sha256"]:
            raise ValueError(f"SHA-256 mismatch for {path}")


def sample_metadata(path: Path = METADATA) -> dict[str, dict[str, str]]:
    result = {}
    for sample in parse_soft_samples(path):
        characteristics = sample["characteristics"]
        result[sample["title"]] = {
            "sample_accession": sample["sample_accession"],
            "stage": sample.get("source", ""),
            "tissue": characteristics.get("tissue", ""),
            "genotype": characteristics.get("genotype", ""),
            "sex": characteristics.get("treatment", ""),
        }
    return result


def build_contrasts(
    headers: list[str],
    metadata: dict[str, dict[str, str]],
    sex: str = "",
) -> list[dict]:
    contrasts = []
    for stage in STAGES:
        for tissue in TISSUES:
            wt = [
                index for index, title in enumerate(headers)
                if metadata[title]["stage"] == stage
                and metadata[title]["tissue"] == tissue
                and metadata[title]["genotype"] == "WT"
                and (not sex or metadata[title]["sex"] == sex)
            ]
            for knockout in KNOCKOUTS:
                ko = [
                    index for index, title in enumerate(headers)
                    if metadata[title]["stage"] == stage
                    and metadata[title]["tissue"] == tissue
                    and metadata[title]["genotype"] == knockout
                    and (not sex or metadata[title]["sex"] == sex)
                ]
                if not wt or not ko:
                    continue
                contrasts.append({
                    "contrast": (
                        f"{stage}_{tissue}_{knockout}_vs_WT"
                        + (f"_{sex}" if sex else "")
                    ),
                    "stage": stage,
                    "tissue": tissue,
                    "knockout": knockout,
                    "sex": sex or "pooled",
                    "WT_indexes": wt,
                    "KO_indexes": ko,
                    "deltas": [],
                    "top_hypomethylated": [],
                })
    return contrasts


def mean_present(values: list[float | None], indexes: list[int]) -> float | None:
    present = [values[index] for index in indexes if values[index] is not None]
    return sum(present) / len(present) if present else None


def parse_beta(value: str) -> float | None:
    if value in {"", "NA", "NaN", "nan"}:
        return None
    parsed = float(value)
    if not 0.0 <= parsed <= 1.0:
        raise ValueError(f"beta value outside [0,1]: {value}")
    return parsed


def analyze(
    matrix: Path = MATRIX,
    metadata_path: Path = METADATA,
) -> dict:
    metadata = sample_metadata(metadata_path)
    sample_sums: defaultdict[str, float] = defaultdict(float)
    sample_counts: defaultdict[str, int] = defaultdict(int)

    with gzip.open(matrix, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        headers = next(reader)
        if set(headers) != set(metadata):
            missing_metadata = sorted(set(headers) - set(metadata))
            missing_matrix = sorted(set(metadata) - set(headers))
            raise ValueError(
                f"matrix/metadata title mismatch: {missing_metadata=}, "
                f"{missing_matrix=}"
            )
        contrasts = build_contrasts(headers, metadata)
        sex_contrasts = (
            build_contrasts(headers, metadata, sex="female")
            + build_contrasts(headers, metadata, sex="male")
        )
        probe_count = 0
        for row in reader:
            if len(row) != len(headers) + 1:
                raise ValueError(
                    f"expected {len(headers) + 1} columns, found {len(row)}"
                )
            probe_id = row[0]
            values = [parse_beta(value) for value in row[1:]]
            probe_count += 1
            for index, value in enumerate(values):
                if value is not None:
                    sample_sums[headers[index]] += value
                    sample_counts[headers[index]] += 1
            for contrast in contrasts + sex_contrasts:
                wt = mean_present(values, contrast["WT_indexes"])
                ko = mean_present(values, contrast["KO_indexes"])
                if wt is None or ko is None:
                    continue
                delta = ko - wt
                contrast["deltas"].append(delta)
                contrast["top_hypomethylated"].append((delta, probe_id, wt, ko))

    sample_rows = []
    for title in headers:
        row = metadata[title]
        sample_rows.append({
            "sample_title": title,
            **row,
            "probes_present": sample_counts[title],
            "mean_beta": sample_sums[title] / sample_counts[title],
        })

    contrast_rows = []
    sex_contrast_rows = []
    top_rows = []
    for contrast in contrasts + sex_contrasts:
        deltas = contrast.pop("deltas")
        ranked = contrast.pop("top_hypomethylated")
        contrast_row = {
            "contrast": contrast["contrast"],
            "stage": contrast["stage"],
            "tissue": contrast["tissue"],
            "knockout": contrast["knockout"],
            "sex": contrast["sex"],
            "WT_samples": len(contrast["WT_indexes"]),
            "KO_samples": len(contrast["KO_indexes"]),
            "probes_compared": len(deltas),
            "mean_KO_minus_WT": statistics.fmean(deltas),
            "median_KO_minus_WT": statistics.median(deltas),
            "fraction_hypomethylated_ge_0.10": (
                sum(delta <= -EFFECT_THRESHOLD for delta in deltas) / len(deltas)
            ),
            "fraction_hypermethylated_ge_0.10": (
                sum(delta >= EFFECT_THRESHOLD for delta in deltas) / len(deltas)
            ),
        }
        if contrast["sex"] == "pooled":
            contrast_rows.append(contrast_row)
            for delta, probe_id, wt, ko in sorted(ranked)[:200]:
                top_rows.append({
                    "contrast": contrast["contrast"],
                    "probe_id": probe_id,
                    "WT_mean_beta": wt,
                    "KO_mean_beta": ko,
                    "KO_minus_WT": delta,
                })
        else:
            sex_contrast_rows.append(contrast_row)

    return {
        "probe_count": probe_count,
        "sample_rows": sample_rows,
        "contrast_rows": contrast_rows,
        "sex_contrast_rows": sex_contrast_rows,
        "top_rows": top_rows,
    }


def write_results(result: dict, output_dir: Path = RESULTS) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = (
        ("sample_global_beta.tsv", result["sample_rows"]),
        ("developmental_contrasts.tsv", result["contrast_rows"]),
        ("developmental_contrasts_by_sex.tsv", result["sex_contrast_rows"]),
        ("top_hypomethylated_probes.tsv", result["top_rows"]),
    )
    for filename, rows in outputs:
        with (output_dir / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, delimiter="\t", fieldnames=list(rows[0])
            )
            writer.writeheader()
            writer.writerows(rows)

    brain = {
        (row["stage"], row["knockout"]): row
        for row in result["contrast_rows"]
        if row["tissue"] == "brain"
    }
    summary = {
        "study": "GSE295720",
        "assay": "Illumina MM285 methylation array",
        "probes": result["probe_count"],
        "samples": len(result["sample_rows"]),
        "brain_stage_isoform_effects": {
            f"{stage}_{knockout}": brain[(stage, knockout)]["mean_KO_minus_WT"]
            for stage in STAGES
            for knockout in KNOCKOUTS
        },
        "developmental_division_of_labor": {
            "Dnmt3a1_effect_change_PD21_minus_E15.5": (
                brain[("PD21", "Dnmt3a1-/-")]["mean_KO_minus_WT"]
                - brain[("E15.5", "Dnmt3a1-/-")]["mean_KO_minus_WT"]
            ),
            "Dnmt3a2_effect_change_PD21_minus_E15.5": (
                brain[("PD21", "Dnmt3a2-/-")]["mean_KO_minus_WT"]
                - brain[("E15.5", "Dnmt3a2-/-")]["mean_KO_minus_WT"]
            ),
        },
        "inference_limit": (
            "Probe-level effect distributions are descriptive. Formal biological "
            "inference must use animals as replicates and account for sex, array "
            "quality, genomic annotation, and multiple testing."
        ),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main() -> None:
    validate_sources()
    result = analyze()
    write_results(result)
    print(json.dumps({
        "probes": result["probe_count"],
        "samples": len(result["sample_rows"]),
        "brain_contrasts": [
            row for row in result["contrast_rows"] if row["tissue"] == "brain"
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
