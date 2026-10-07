#!/usr/bin/env python3
"""Validate brain phenotypes before comparative association analyses."""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/traits/brain_phenotypes.tsv"
TAXA = ROOT / "results/01_curated/taxon_metadata.tsv"
RESULTS = ROOT / "results/05_brain_integration"
MIN_MAMMAL_SPECIES_FOR_PILOT_ASSOCIATION = 10
REQUIRED = {
    "species", "phenotype", "value", "unit", "brain_region", "cell_type",
    "developmental_stage", "assay", "source_accession", "citation_doi",
}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def validate(
    rows: list[dict[str, str]], taxa: list[dict[str, str]]
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, object]]:
    issues: list[dict[str, str]] = []
    names: dict[str, dict[str, str]] = {}
    for taxon in taxa:
        for key in ("species", "safe_id", "ncbi_scientific_name"):
            names[taxon[key]] = taxon

    seen: set[tuple[str, ...]] = set()
    valid_rows: list[tuple[dict[str, str], dict[str, str]]] = []
    for index, row in enumerate(rows, start=2):
        for field in REQUIRED:
            if not row.get(field, "").strip():
                issues.append({
                    "severity": "error", "row": str(index),
                    "code": "missing_required", "message": field,
                })
        taxon = names.get(row.get("species", ""))
        if taxon is None:
            issues.append({
                "severity": "error", "row": str(index),
                "code": "species_not_in_panel",
                "message": row.get("species", ""),
            })
        try:
            value = float(row.get("value", ""))
            if not math.isfinite(value):
                raise ValueError
            if row.get("unit") == "fraction" and not 0 <= value <= 1:
                issues.append({
                    "severity": "error", "row": str(index),
                    "code": "fraction_out_of_range", "message": str(value),
                })
        except ValueError:
            issues.append({
                "severity": "error", "row": str(index),
                "code": "invalid_numeric_value",
                "message": row.get("value", ""),
            })

        duplicate_key = tuple(
            row.get(field, "") for field in
            ("species", "phenotype", "brain_region", "developmental_stage",
             "sample_id", "source_accession")
        )
        if duplicate_key in seen:
            issues.append({
                "severity": "error", "row": str(index),
                "code": "duplicate_measurement",
                "message": "|".join(duplicate_key),
            })
        seen.add(duplicate_key)
        if taxon is not None:
            valid_rows.append((row, taxon))

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    taxon_by_species: dict[str, dict[str, str]] = {}
    for row, taxon in valid_rows:
        scientific_name = taxon["species"]
        grouped[scientific_name].append(row)
        taxon_by_species[scientific_name] = taxon

    overlap = []
    for species in sorted(grouped):
        taxon = taxon_by_species[species]
        records = grouped[species]
        overlap.append({
            "species": species,
            "safe_id": taxon["safe_id"],
            "is_mammal": taxon["is_mammal"],
            "phenotype_records": str(len(records)),
            "phenotypes": ";".join(sorted({row["phenotype"] for row in records})),
            "source_accessions": ";".join(
                sorted({row["source_accession"] for row in records})
            ),
        })

    error_count = sum(issue["severity"] == "error" for issue in issues)
    mammal_species = {
        taxon["species"] for _, taxon in valid_rows if taxon["is_mammal"] == "true"
    }
    mammals_by_phenotype: dict[str, set[str]] = defaultdict(set)
    for row, taxon in valid_rows:
        if taxon["is_mammal"] == "true":
            mammals_by_phenotype[row["phenotype"]].add(taxon["species"])
    mammal_counts_by_phenotype = {
        phenotype: len(species)
        for phenotype, species in sorted(mammals_by_phenotype.items())
    }
    ready_phenotypes = sorted(
        phenotype for phenotype, count in mammal_counts_by_phenotype.items()
        if count >= MIN_MAMMAL_SPECIES_FOR_PILOT_ASSOCIATION
    )
    phenotype_counts = Counter(row["phenotype"] for row, _ in valid_rows)
    summary: dict[str, object] = {
        "records": len(rows),
        "matched_panel_records": len(valid_rows),
        "matched_panel_species": len(grouped),
        "matched_mammal_species": len(mammal_species),
        "phenotype_counts": dict(sorted(phenotype_counts.items())),
        "mammal_species_by_phenotype": mammal_counts_by_phenotype,
        "errors": error_count,
        "warnings": sum(issue["severity"] == "warning" for issue in issues),
        "minimum_mammal_species_for_pilot_association":
            MIN_MAMMAL_SPECIES_FOR_PILOT_ASSOCIATION,
        "pilot_ready_phenotypes": ready_phenotypes,
        "pilot_phylogenetic_association_ready": (
            error_count == 0 and bool(ready_phenotypes)
        ),
    }
    if not issues:
        issues.append({
            "severity": "info", "row": "",
            "code": "validation_passed",
            "message": f"{len(rows)} records passed schema and value checks",
        })
    return issues, overlap, summary


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    fields, rows = read_tsv(INPUT)
    missing_columns = sorted(REQUIRED - set(fields))
    if missing_columns:
        raise SystemExit("Missing phenotype columns: " + ", ".join(missing_columns))
    _, taxa = read_tsv(TAXA)
    issues, overlap, summary = validate(rows, taxa)
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_tsv(
        RESULTS / "phenotype_validation.tsv", issues,
        ["severity", "row", "code", "message"],
    )
    write_tsv(
        RESULTS / "phenotype_species_overlap.tsv", overlap,
        ["species", "safe_id", "is_mammal", "phenotype_records",
         "phenotypes", "source_accessions"],
    )
    (RESULTS / "phenotype_validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))
    if summary["errors"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
