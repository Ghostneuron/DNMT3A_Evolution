#!/usr/bin/env python3
"""Summarize tissue context of DNMT3A internal-promoter FANTOM5 counts.

Counts are used for within-library detection only. They are not compared as
expression magnitudes across libraries because this module does not normalize
library depth.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/fantom5_cage"
RESULTS = ROOT / "results/06_isoform_evolution"
DIRECT = RESULTS / "direct_tss_evidence.tsv"
NEARBY = RESULTS / "fantom5_nearby_cage_peaks.tsv"
MATRICES = {
    "Homo sapiens": RAW / "hg38_fair_new_CAGE_peaks_phase1and2_counts.osc.txt.gz",
    "Mus musculus": RAW / "mm10_fair_new_CAGE_peaks_phase1and2_counts.osc.txt.gz",
}
BRAIN_PATTERNS = {
    "brain": re.compile(
        r"\b(brain|cerebell\w*|cerebral|cortex|cortical|hippocamp\w*|"
        r"amygdala|thalam\w*|striat\w*|substantia nigra|medulla oblongata|"
        r"spinal cord|astrocyte\w*|oligodendro\w*|microglia\w*|"
        r"neural stem|neural progenitor|neuron\w*)\b",
        re.I,
    ),
    "non_neural_cortex": re.compile(r"\b(renal|kidney|adrenal)\b", re.I),
    "developmental": re.compile(
        r"\b(fetal|foetal|embry\w*|newborn|neonat\w*|postnatal|"
        r"day\s*\d+|week\s*\d+|month\s*\d+|differentiation|"
        r"progenitor|stem cell)\b",
        re.I,
    ),
}
SAMPLE_ACCESSION = re.compile(r"\.(CNhs\d+)\.")


def decode_sample(column: str) -> tuple[str, str]:
    accession_match = SAMPLE_ACCESSION.search(column)
    accession = accession_match.group(1) if accession_match else ""
    encoded = column.removeprefix("counts.")
    if accession_match:
        encoded = encoded[:encoded.find(f".{accession}.")]
    return accession, unquote(encoded)


def tissue_flags(description: str) -> tuple[bool, bool]:
    brain = bool(BRAIN_PATTERNS["brain"].search(description))
    if BRAIN_PATTERNS["non_neural_cortex"].search(description):
        brain = False
    return (
        brain,
        bool(BRAIN_PATTERNS["developmental"].search(description)),
    )


def read_target_counts(path: Path, target_ids: set[str]):
    columns = None
    target_rows = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("00Annotation\t"):
                columns = line.rstrip().split("\t")[1:]
            elif not line.startswith("#") and columns is not None:
                fields = line.rstrip().split("\t")
                if fields[0] in target_ids:
                    target_rows[fields[0]] = [int(value) for value in fields[1:]]
                    if len(target_rows) == len(target_ids):
                        break
    if columns is None:
        raise ValueError(f"No OSC column header in {path}")
    missing = target_ids - set(target_rows)
    if missing:
        raise ValueError(f"Missing target peaks in {path}: {sorted(missing)}")
    if any(len(values) != len(columns) for values in target_rows.values()):
        raise ValueError(f"Count/header length mismatch in {path}")
    return columns, target_rows


def main() -> None:
    with DIRECT.open(newline="") as handle:
        direct = {
            row["species"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["direct_tss_source"] == "FANTOM5 CAGE"
        }
    with NEARBY.open(newline="") as handle:
        nearby = list(csv.DictReader(handle, delimiter="\t"))

    sample_rows = []
    summaries = []
    for species, matrix in MATRICES.items():
        local_peaks = {
            row["peak_id"] for row in nearby if row["species"] == species
        }
        closest_peak = direct[species]["closest_peak_id"]
        columns, counts = read_target_counts(matrix, local_peaks)
        for index, column in enumerate(columns):
            accession, description = decode_sample(column)
            brain, developmental = tissue_flags(description)
            closest_count = counts[closest_peak][index]
            cluster_count = sum(values[index] for values in counts.values())
            sample_rows.append({
                "species": species,
                "sample_accession": accession,
                "sample_description": description,
                "brain_CNS_keyword_match": str(brain).lower(),
                "developmental_keyword_match": str(developmental).lower(),
                "closest_CAGE_peak_id": closest_peak,
                "closest_peak_raw_count": closest_count,
                "local_500bp_peak_count": len(local_peaks),
                "local_500bp_cluster_raw_count": cluster_count,
                "closest_peak_detected": str(closest_count > 0).lower(),
                "local_cluster_detected": str(cluster_count > 0).lower(),
                "quantification_limit": (
                    "Raw counts used for detection; magnitudes are not "
                    "compared across unnormalized libraries"
                ),
            })

        species_rows = [row for row in sample_rows if row["species"] == species]
        brain_rows = [
            row for row in species_rows if row["brain_CNS_keyword_match"] == "true"
        ]
        nonbrain_rows = [
            row for row in species_rows if row["brain_CNS_keyword_match"] == "false"
        ]
        developmental_brain = [
            row for row in brain_rows
            if row["developmental_keyword_match"] == "true"
        ]
        summaries.append({
            "species": species,
            "all_samples": len(species_rows),
            "brain_CNS_keyword_samples": len(brain_rows),
            "developmental_brain_keyword_samples": len(developmental_brain),
            "closest_peak_detected_all_samples": detected(species_rows, "closest_peak_detected"),
            "closest_peak_detected_brain_samples": detected(brain_rows, "closest_peak_detected"),
            "closest_peak_detected_developmental_brain_samples": detected(
                developmental_brain, "closest_peak_detected"
            ),
            "closest_peak_detected_nonbrain_samples": detected(nonbrain_rows, "closest_peak_detected"),
            "local_cluster_detected_brain_samples": detected(brain_rows, "local_cluster_detected"),
            "local_cluster_detected_developmental_brain_samples": detected(
                developmental_brain, "local_cluster_detected"
            ),
            "brain_classification_method": "transparent sample-description keyword screen",
        })

    with (RESULTS / "fantom5_internal_promoter_sample_counts.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(sample_rows[0])
        )
        writer.writeheader()
        writer.writerows(sample_rows)

    summary = {
        "species": summaries,
        "brain_specificity_inference": False,
        "interpretation": (
            "Detection in keyword-classified brain/CNS samples establishes "
            "brain-context activity. It does not establish enrichment or "
            "developmental specificity without normalized expression and "
            "a curated sample ontology."
        ),
    }
    (RESULTS / "fantom5_tissue_context_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


def detected(rows: list[dict[str, str]], field: str) -> dict[str, float | int | None]:
    count = sum(row[field] == "true" for row in rows)
    return {
        "detected": count,
        "total": len(rows),
        "fraction": count / len(rows) if rows else None,
    }


if __name__ == "__main__":
    main()
