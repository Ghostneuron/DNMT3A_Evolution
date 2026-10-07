#!/usr/bin/env python3
"""Classify annotated DNMT3A products relative to each species' full form."""

from __future__ import annotations

import csv
import difflib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from dnmt3a_pipeline import read_fasta


ROOT = Path(__file__).resolve().parents[1]
PROTEINS = ROOT / "data/raw/ncbi_dataset/data/protein.faa"
REPRESENTATIVES = ROOT / "results/01_curated/representatives.tsv"
TAXA = ROOT / "results/01_curated/taxon_metadata.tsv"
RESULTS = ROOT / "results/06_isoform_evolution"

ORGANISM = re.compile(r"\[organism=([^\]]+)\]")
ISOFORM = re.compile(r"\[isoform=([^\]]+)\]")


def longest_shared_block(full: str, alternative: str) -> tuple[int, int, int]:
    match = difflib.SequenceMatcher(
        None, full, alternative, autojunk=False
    ).find_longest_match()
    return match.a, match.b, match.size


def classify_product(
    full: str, alternative: str, selected: bool = False
) -> tuple[str, str]:
    if selected:
        return "selected_full_length", "selected coding-evolution representative"
    full_start, alternative_start, shared = longest_shared_block(full, alternative)
    alternative_coverage = shared / len(alternative)
    full_fraction = full_start / len(full)
    shared_fraction = shared / len(full)

    downstream_core = (
        0.18 <= full_fraction <= 0.28
        and shared_fraction >= 0.70
        and alternative_coverage >= 0.90
    )
    if downstream_core and 8 <= alternative_start <= 45:
        return (
            "DNMT3A2_like_downstream_start",
            "unique short leader followed by the conserved downstream core",
        )
    if downstream_core:
        if alternative_start < 8:
            return (
                "leaderless_downstream_core_product",
                "exact or near-exact downstream core without an A2-like leader",
            )
        return (
            "extended_leader_downstream_core_product",
            "downstream shared core preceded by a leader longer than the A2 rule",
        )
    if len(alternative) >= 0.85 * len(full) and alternative_coverage >= 0.85:
        return "DNMT3A1_like_alternative", "near-full-length alternative product"
    if len(alternative) < 300 and full_start < 30:
        return "N_terminal_fragment", "short product dominated by the DNMT3A1 N terminus"
    return "other_alternative_product", "does not meet operational A1/A2 rules"


def evidence_tier(accession: str) -> str:
    if accession.startswith("NP_"):
        return "curated_RefSeq"
    if accession.startswith("XP_"):
        return "predicted_RefSeq"
    return "other"


def mammal_group(lineage: str) -> str:
    if "Prototheria" in lineage:
        return "Prototheria"
    if "Metatheria" in lineage:
        return "Metatheria"
    return "Eutheria"


def main() -> None:
    with REPRESENTATIVES.open(newline="") as handle:
        representatives = list(csv.DictReader(handle, delimiter="\t"))
    selected_by_species = {
        row["species"]: row["protein_accession"] for row in representatives
    }
    with TAXA.open(newline="") as handle:
        taxa = {
            row["species"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["is_mammal"] == "true"
        }

    proteins_by_species = defaultdict(list)
    for record in read_fasta(PROTEINS):
        match = ORGANISM.search(record.description)
        if match and match.group(1) in taxa:
            proteins_by_species[match.group(1)].append(record)

    rows = []
    for species in sorted(taxa):
        selected_accession = selected_by_species[species]
        records = proteins_by_species[species]
        by_accession = {record.identifier: record for record in records}
        full_record = by_accession[selected_accession]
        for record in records:
            full_start, product_start, shared = longest_shared_block(
                full_record.sequence, record.sequence
            )
            product_class, reason = classify_product(
                full_record.sequence,
                record.sequence,
                selected=record.identifier == selected_accession,
            )
            isoform_match = ISOFORM.search(record.description)
            rows.append({
                "species": species,
                "safe_id": taxa[species]["safe_id"],
                "order": taxa[species]["order"],
                "mammal_group": mammal_group(taxa[species]["lineage"]),
                "protein_accession": record.identifier,
                "annotated_isoform": isoform_match.group(1) if isoform_match else "",
                "evidence_tier": evidence_tier(record.identifier),
                "protein_length_aa": len(record.sequence),
                "selected_full_length_aa": len(full_record.sequence),
                "shared_block_full_start_0based": full_start,
                "shared_block_product_start_0based": product_start,
                "shared_block_length_aa": shared,
                "shared_block_product_coverage": shared / len(record.sequence),
                "product_class": product_class,
                "classification_reason": reason,
            })

    fields = list(rows[0])
    RESULTS.mkdir(parents=True, exist_ok=True)
    with (RESULTS / "isoform_product_classification.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    by_species = defaultdict(list)
    for row in rows:
        by_species[row["species"]].append(row)
    species_rows = []
    for species in sorted(by_species):
        products = by_species[species]
        candidates = [
            row for row in products
            if row["product_class"] == "DNMT3A2_like_downstream_start"
        ]
        curated = [
            row for row in candidates if row["evidence_tier"] == "curated_RefSeq"
        ]
        if curated:
            status = "curated_A2_like_annotation_present"
        elif candidates:
            status = "predicted_A2_like_annotation_present"
        else:
            status = "A2_like_not_detected_in_NCBI_package"
        leaderless = [
            row for row in products
            if row["product_class"] == "leaderless_downstream_core_product"
        ]
        extended = [
            row for row in products
            if row["product_class"] == "extended_leader_downstream_core_product"
        ]
        species_rows.append({
            "species": species,
            "safe_id": taxa[species]["safe_id"],
            "order": taxa[species]["order"],
            "mammal_group": mammal_group(taxa[species]["lineage"]),
            "annotated_protein_products": len(products),
            "unique_protein_sequences": len({
                next(
                    record.sequence for record in proteins_by_species[species]
                    if record.identifier == row["protein_accession"]
                )
                for row in products
            }),
            "A2_like_products": len(candidates),
            "curated_A2_like_products": len(curated),
            "leaderless_downstream_core_products": len(leaderless),
            "extended_leader_downstream_core_products": len(extended),
            "A2_like_accessions": ";".join(
                row["protein_accession"] for row in candidates
            ),
            "annotation_status": status,
            "broader_downstream_core_status": (
                "strict_A2_like_present" if candidates else
                "leaderless_core_present" if leaderless else
                "extended_leader_core_present" if extended else
                "no_downstream_core_product_detected"
            ),
            "absence_interpretation": (
                "" if candidates else
                "annotation non-detection; not evidence of biological absence"
            ),
        })
    species_fields = list(species_rows[0])
    with (RESULTS / "isoform_species_summary.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=species_fields
        )
        writer.writeheader()
        writer.writerows(species_rows)

    class_counts = Counter(row["product_class"] for row in rows)
    status_counts = Counter(row["annotation_status"] for row in species_rows)
    human = [
        row for row in rows
        if row["species"] == "Homo sapiens"
        and row["protein_accession"] in {"NP_072046.2", "NP_715640.2"}
    ]
    summary = {
        "mammal_species": len(species_rows),
        "annotated_protein_products": len(rows),
        "product_class_counts": dict(sorted(class_counts.items())),
        "species_annotation_status_counts": dict(sorted(status_counts.items())),
        "human_calibration": {
            row["protein_accession"]: row["product_class"] for row in human
        },
        "inference_limit": (
            "NCBI annotation presence supports an expressed/predicted product; "
            "non-detection does not establish isoform loss"
        ),
    }
    (RESULTS / "isoform_classification_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
