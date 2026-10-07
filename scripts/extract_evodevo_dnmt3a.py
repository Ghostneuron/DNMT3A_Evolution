#!/usr/bin/env python3
"""Extract a DNMT3A row from an EvoDevo CPM matrix and make brain traits."""

from __future__ import annotations

import argparse
import csv
import re
import shlex
from collections import defaultdict
from pathlib import Path

try:
    from .import_brainzoo_mch import FIELDS, ROOT
except ImportError:  # Direct script execution.
    from import_brainzoo_mch import FIELDS, ROOT


DOI = "10.1038/s41586-019-1338-5"


def extract_row(matrix: Path, gene_id: str) -> tuple[list[str], list[str]]:
    with matrix.open() as handle:
        samples = shlex.split(handle.readline())
        for line in handle:
            if line.startswith(f'"{gene_id}"') or line.startswith(gene_id + " "):
                values = shlex.split(line)
                if len(values) != len(samples) + 1:
                    raise ValueError(
                        f"{gene_id}: {len(values) - 1} values for "
                        f"{len(samples)} samples"
                    )
                return samples, values[1:]
    raise ValueError(f"{gene_id} not found in {matrix}")


def load_sdrf_age_map(path: Path) -> dict[str, tuple[str, str]]:
    """Return unambiguous age -> (developmental stage, unit) mappings."""
    observations: dict[str, set[tuple[str, str]]] = defaultdict(set)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            age = row.get("Characteristics[age]", "").strip()
            stage = row.get("Characteristics[developmental stage]", "").strip()
            unit = row.get("Unit[time unit]", "").strip()
            if age and stage and unit:
                observations[age].add((stage, unit))
    return {
        age: next(iter(values))
        for age, values in observations.items() if len(values) == 1
    }


def stage_fields(
    stage: str,
    age_map: dict[str, tuple[str, str]] | None = None,
    postnatal_offset_days: float | None = None,
) -> tuple[str, str, str]:
    match = re.fullmatch(r"(\d+(?:\.\d+)?)wpc", stage)
    if match:
        return "embryo", match.group(1), "week"
    match = re.fullmatch(r"(\d+(?:\.\d+)?)wpb", stage)
    if match:
        return "postnatal", match.group(1), "week"
    match = re.fullmatch(r"(\d+(?:\.\d+)?)mpb", stage)
    if match:
        return "postnatal", match.group(1), "month"
    match = re.fullmatch(r"[eE](\d+(?:\.\d+)?)", stage)
    if match:
        return "embryo", match.group(1), "day"
    match = re.fullmatch(r"[pP](\d+(?:\.\d+)?)", stage)
    if match:
        return "postnatal", match.group(1), "day"
    if age_map and postnatal_offset_days is not None:
        direct = age_map.get(stage)
        if direct and direct[0] == "embryo":
            developmental_stage, unit = direct
            return developmental_stage, stage, unit
        try:
            adjusted = float(stage) - postnatal_offset_days
            adjusted_label = f"{adjusted:g}"
        except ValueError:
            adjusted_label = ""
        if adjusted_label in age_map:
            developmental_stage, unit = age_map[adjusted_label]
            return developmental_stage, adjusted_label, unit
    if age_map and stage in age_map:
        developmental_stage, unit = age_map[stage]
        return developmental_stage, stage, unit
    return stage, "", ""


def brain_rows(
    samples: list[str],
    values: list[str],
    species: str,
    accession: str,
    gene_id: str,
    age_map: dict[str, tuple[str, str]] | None = None,
    postnatal_offset_days: float | None = None,
) -> list[dict[str, str]]:
    rows = []
    for sample, value in zip(samples, values):
        parts = sample.split(".")
        if len(parts) < 3 or parts[0] not in {"Brain", "Cerebellum"}:
            continue
        region, stage = parts[0], ".".join(parts[1:-1])
        developmental_stage, age_value, age_unit = stage_fields(
            stage, age_map, postnatal_offset_days
        )
        rows.append({
            "species": species,
            "phenotype": "brain_total_DNMT3A_expression",
            "value": value,
            "unit": "CPM",
            "brain_region": (
                "brain/cerebrum series" if region == "Brain" else "cerebellum"
            ),
            "cell_type": "bulk tissue",
            "developmental_stage": developmental_stage,
            "age_value": age_value,
            "age_unit": age_unit,
            "sex": "",
            "assay": "RNA-seq",
            "sample_id": sample,
            "source_accession": accession,
            "source_table": (
                f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}"
            ),
            "citation_doi": DOI,
            "normalization": "TMM-normalized counts per million",
            "derivation": f"processed CPM row {gene_id}",
            "quality_notes": (
                "Gene-level total DNMT3A; not DNMT3A1/DNMT3A2 isoform-specific. "
                "Brain series changes anatomical definition across development. "
                f"Processed matrix stage label: {stage}."
            ),
        })
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--sdrf", type=Path)
    parser.add_argument(
        "--postnatal-offset-days",
        type=float,
        help=(
            "Subtract this gestation offset from numeric processed stage labels "
            "before matching postnatal SDRF ages"
        ),
    )
    parser.add_argument("--raw-subset", type=Path, required=True)
    parser.add_argument("--trait-output", type=Path, required=True)
    args = parser.parse_args()

    samples, values = extract_row(args.matrix, args.gene_id)
    write_tsv(
        args.raw_subset,
        [
            {"gene_id": args.gene_id, "sample_id": sample, "CPM": value}
            for sample, value in zip(samples, values)
        ],
        ["gene_id", "sample_id", "CPM"],
    )
    age_map = load_sdrf_age_map(args.sdrf) if args.sdrf else None
    rows = brain_rows(
        samples, values, args.species, args.accession, args.gene_id, age_map,
        args.postnatal_offset_days,
    )
    write_tsv(args.trait_output, rows, FIELDS)
    print(
        f"Extracted {len(samples)} total samples and {len(rows)} brain records "
        f"for {args.species}"
    )


if __name__ == "__main__":
    main()
