#!/usr/bin/env python3
"""Import BrainZoo global brain non-CG methylation summaries."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/raw/brain/BrainZoo_mCH_stats.tsv"
TAXA = ROOT / "results/01_curated/taxon_metadata.tsv"
OUTPUT = ROOT / "data/traits/sources/brainzoo_mch.tsv"
RESULTS = ROOT / "results/05_brain_integration"
DOI = "10.1038/s41559-020-01371-2"
SOURCE_ID = "BrainZoo:mCH_stats.tsv"
SOURCE_URL = "https://github.com/AlexdeMendoza/BrainZoo/blob/master/mCH_stats.tsv"

SPECIES_MAP = {
    "Shark": "Callorhinchus milii",
    "Chicken": "Gallus gallus",
    "Lamprey": "Lethenteron camtschaticum",
    "Mouse": "Mus musculus",
    "Opossum": "Monodelphis domestica",
    "Platypus": "Ornithorhynchus anatinus",
    "Octopus": "Octopus bimaculoides",
    "Xenopus": "Xenopus laevis",
    "Honeybee": "Apis mellifera",
    "Human": "Homo sapiens",
    "Parus": "Parus major",
    "Zebrafish_forebrain": "Danio rerio",
    "Lancelet_neural": "Branchiostoma lanceolatum",
}

BRAIN_REGION = {
    "Human": "middle frontal gyrus",
    "Mouse": "frontal cortex",
    "Zebrafish_forebrain": "forebrain",
    "Octopus": "supraesophageal brain",
    "Lancelet_neural": "neural tube",
}

FIELDS = [
    "species", "phenotype", "value", "unit", "brain_region", "cell_type",
    "developmental_stage", "age_value", "age_unit", "sex", "assay",
    "sample_id", "source_accession", "source_table", "citation_doi",
    "normalization", "derivation", "quality_notes",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_brainzoo(path: Path) -> list[dict[str, str]]:
    """Parse the source's R-style row-name column without shifting fields."""
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = []
        for line_number, line in enumerate(handle, start=2):
            values = line.rstrip("\n").split("\t")
            if len(values) == len(header) + 1:
                values = values[1:]
            if len(values) != len(header):
                raise ValueError(
                    f"{path}:{line_number}: {len(values)} fields; "
                    f"expected {len(header)}"
                )
            rows.append(dict(zip(header, values)))
    return rows


def phenotype_rows(
    source_rows: list[dict[str, str]], panel_species: set[str]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    output: list[dict[str, str]] = []
    audit: list[dict[str, str]] = []
    for row in source_rows:
        label = row["species"]
        scientific_name = SPECIES_MAP.get(label, "")
        status = "matched_panel" if scientific_name in panel_species else "not_in_panel"
        if not scientific_name:
            status = "unmapped_label"
        audit.append({
            "source_label": label,
            "scientific_name": scientific_name,
            "sample_id": row["file"],
            "panel_status": status,
        })
        if status != "matched_panel":
            continue

        measurements = {
            "global_brain_mCH_fraction": (
                int(row["total_mCH"]) / int(row["total_CHcov"]),
                "total_mCH / total_CHcov",
            ),
            "global_brain_mCA_fraction": (float(row["mCA"]), "source mCA"),
            "global_brain_mCT_fraction": (float(row["mCT"]), "source mCT"),
            "global_brain_mCC_fraction": (float(row["mCC"]), "source mCC"),
        }
        for phenotype, (value, derivation) in measurements.items():
            output.append({
                "species": scientific_name,
                "phenotype": phenotype,
                "value": f"{value:.15g}",
                "unit": "fraction",
                "brain_region": BRAIN_REGION.get(label, "brain"),
                "cell_type": "bulk tissue",
                "developmental_stage": "adult",
                "age_value": "25" if label == "Human" else ("6" if label == "Mouse" else ""),
                "age_unit": "years" if label == "Human" else ("weeks" if label == "Mouse" else ""),
                "sex": "",
                "assay": "WGBS",
                "sample_id": row["file"],
                "source_accession": SOURCE_ID,
                "source_table": SOURCE_URL,
                "citation_doi": DOI,
                "normalization": "coverage-weighted global methylated fraction",
                "derivation": derivation,
                "quality_notes": (
                    "Compact published analysis summary; generated and "
                    "reanalyzed WGBS samples have heterogeneous source studies"
                ),
            })
    return output, audit


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    panel_species = {row["species"] for row in read_tsv(TAXA)}
    source_rows = parse_brainzoo(SOURCE)
    rows, audit = phenotype_rows(source_rows, panel_species)
    write_tsv(OUTPUT, rows, FIELDS)
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_tsv(
        RESULTS / "brainzoo_import_audit.tsv",
        audit,
        ["source_label", "scientific_name", "sample_id", "panel_status"],
    )
    matched_species = sorted({row["species"] for row in rows})
    summary = {
        "source_rows": len(source_rows),
        "matched_panel_species": len(matched_species),
        "matched_species": matched_species,
        "phenotype_records": len(rows),
        "measurements_per_species": 4,
    }
    (RESULTS / "brainzoo_import_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
