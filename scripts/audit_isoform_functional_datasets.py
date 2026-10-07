#!/usr/bin/env python3
"""Inventory public functional-genomics datasets that distinguish DNMT3A isoforms."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/isoform_functional_genomics"
RESULTS = ROOT / "results/07_isoform_function"
STUDIES = {
    "GSE164265": {
        "citation": "Gu et al. 2022, Nature Genetics",
        "primary_use": (
            "P21 cortex and sorted cortical-neuron methylome/transcriptome; "
            "isoform-specific knockout and DNMT3A1-tail mechanism"
        ),
        "comparison_tier": "direct_matched_brain_isoform_comparison",
    },
    "GSE295720": {
        "citation": "Liu et al. 2025, Communications Biology",
        "primary_use": (
            "E15.5 versus P21 brain/liver methylation-array comparison of "
            "Dnmt3a1 and Dnmt3a2 knockouts"
        ),
        "comparison_tier": "direct_matched_developmental_isoform_comparison",
    },
    "GSE96529": {
        "citation": "Manzo et al. 2017, EMBO Journal",
        "primary_use": (
            "Controlled DNMT3A1/DNMT3A2 rescue, binding, and 5hmC assays "
            "in mouse embryonic stem cells"
        ),
        "comparison_tier": "direct_matched_mechanistic_isoform_comparison",
    },
}


def parse_soft_samples(path: Path) -> list[dict]:
    samples: list[dict] = []
    current: dict | None = None
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                current = {
                    "sample_accession": line.split(" = ", 1)[1],
                    "characteristics": {},
                    "supplementary_files": [],
                }
            elif current is None:
                continue
            elif line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_source_name_ch1 = "):
                current["source"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_library_strategy = "):
                current["library_strategy"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current["characteristics"][key.lower()] = item
            elif line.startswith("!Sample_supplementary_file = "):
                current["supplementary_files"].append(line.split(" = ", 1)[1])
    if current is not None:
        samples.append(current)
    return samples


def flatten_sample(study: str, sample: dict) -> dict[str, str | int]:
    characteristics = sample["characteristics"]
    return {
        "study": study,
        "sample_accession": sample["sample_accession"],
        "title": sample.get("title", ""),
        "source": sample.get("source", ""),
        "tissue": characteristics.get("tissue", ""),
        "genotype": characteristics.get("genotype", ""),
        "developmental_stage": characteristics.get(
            "developmental stage", sample.get("source", "")
        ),
        "treatment_or_sex": characteristics.get("treatment", ""),
        "library_strategy": sample.get("library_strategy", ""),
        "supplementary_file_count": len(sample["supplementary_files"]),
    }


def design_counts(samples: list[dict[str, str | int]]) -> list[dict]:
    keys = (
        "study",
        "source",
        "tissue",
        "developmental_stage",
        "genotype",
        "library_strategy",
    )
    counts = Counter(tuple(row[key] for key in keys) for row in samples)
    return [
        {**dict(zip(keys, values)), "samples": count}
        for values, count in sorted(counts.items())
    ]


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    all_samples = []
    registry = []
    for study, metadata in STUDIES.items():
        path = RAW / study / f"{study}_samples.soft.txt"
        samples = parse_soft_samples(path)
        flattened = [flatten_sample(study, sample) for sample in samples]
        all_samples.extend(flattened)
        registry.append({
            "study": study,
            **metadata,
            "geo_samples": len(samples),
            "metadata_file": str(path.relative_to(ROOT)),
        })

    sample_fields = list(all_samples[0])
    with (RESULTS / "functional_dataset_samples.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=sample_fields)
        writer.writeheader()
        writer.writerows(all_samples)

    counts = design_counts(all_samples)
    with (RESULTS / "functional_dataset_design_counts.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(counts[0])
        )
        writer.writeheader()
        writer.writerows(counts)

    with (RESULTS / "functional_dataset_registry.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(registry[0])
        )
        writer.writeheader()
        writer.writerows(registry)

    summary = {
        "studies": len(registry),
        "samples": len(all_samples),
        "samples_by_study": dict(Counter(row["study"] for row in all_samples)),
        "direct_matched_designs": {
            row["study"]: row["comparison_tier"] for row in registry
        },
        "analysis_priority": [
            "GSE164265 P21 sorted cortical-neuron CpG methylation",
            "GSE295720 E15.5-to-P21 brain methylation-array dynamics",
            "GSE96529 controlled isoform localization and rescue",
        ],
        "cross_study_limit": (
            "Effect sizes must first be estimated within each study. Raw levels "
            "from different assays, tissues, stages, or genome builds are not "
            "directly comparable."
        ),
    }
    (RESULTS / "functional_dataset_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
