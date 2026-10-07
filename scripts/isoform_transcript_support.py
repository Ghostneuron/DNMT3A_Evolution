#!/usr/bin/env python3
"""Link downstream-core protein products to packaged CDS/RNA models."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from dnmt3a_pipeline import header_fields, read_fasta, translate


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/raw/ncbi_dataset/data"
PRODUCTS = ROOT / "results/06_isoform_evolution/isoform_product_classification.tsv"
RESULTS = ROOT / "results/06_isoform_evolution"
CDS_COORDS = re.compile(r"^([^:]+):(\d+)-(\d+)$")
TARGET_CLASSES = {
    "DNMT3A2_like_downstream_start",
    "leaderless_downstream_core_product",
    "extended_leader_downstream_core_product",
}


def accession_tier(accession: str) -> str:
    if accession.startswith("NM_"):
        return "curated_RefSeq_transcript"
    if accession.startswith("XM_"):
        return "predicted_RefSeq_transcript"
    return "other_transcript"


def main() -> None:
    proteins = {
        record.identifier: record
        for record in read_fasta(DATA / "protein.faa")
    }
    rna_lengths = {
        record.identifier: len(record.sequence)
        for record in read_fasta(DATA / "rna.fna")
    }
    cds_by_species_translation = defaultdict(list)
    for record in read_fasta(DATA / "cds.fna"):
        match = CDS_COORDS.match(record.identifier)
        if not match:
            continue
        transcript, start, end = match.groups()
        species = header_fields(record.description).get("organism", "")
        cds_by_species_translation[(species, translate(record.sequence))].append({
            "transcript_accession": transcript,
            "cds_start_1based": int(start),
            "cds_end_1based": int(end),
        })

    with PRODUCTS.open(newline="") as handle:
        candidates = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["product_class"] in TARGET_CLASSES
        ]

    rows = []
    unsupported = []
    for candidate in candidates:
        protein = proteins[candidate["protein_accession"]]
        matches = cds_by_species_translation.get(
            (candidate["species"], protein.sequence), []
        )
        if not matches:
            unsupported.append(candidate["protein_accession"])
        for match in matches:
            transcript = match["transcript_accession"]
            transcript_length = rna_lengths.get(transcript)
            cds_end = match["cds_end_1based"]
            rows.append({
                "species": candidate["species"],
                "mammal_group": candidate["mammal_group"],
                "protein_accession": candidate["protein_accession"],
                "product_class": candidate["product_class"],
                "protein_evidence_tier": candidate["evidence_tier"],
                "transcript_accession": transcript,
                "transcript_evidence_tier": accession_tier(transcript),
                "cds_start_1based": match["cds_start_1based"],
                "cds_end_1based": cds_end,
                "transcript_length_nt": transcript_length or "",
                "five_prime_UTR_length_nt": match["cds_start_1based"] - 1,
                "three_prime_UTR_length_nt": (
                    transcript_length - cds_end if transcript_length else ""
                ),
                "exact_CDS_translation_match": True,
                "evidence_limit": (
                    "Internal NCBI model consistency; XM/XP models are not "
                    "experimental transcription-start evidence"
                ),
            })

    fields = list(rows[0])
    with (RESULTS / "isoform_transcript_model_support.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    non_eutherian = [
        row for row in rows if row["mammal_group"] != "Eutheria"
    ]
    summary = {
        "downstream_core_protein_products": len(candidates),
        "products_with_exact_packaged_CDS_match": len({
            row["protein_accession"] for row in rows
        }),
        "products_without_exact_packaged_CDS_match": len(unsupported),
        "exact_species_matched_protein_CDS_links": len(rows),
        "transcript_model_tiers": dict(sorted(Counter(
            row["transcript_evidence_tier"] for row in rows
        ).items())),
        "non_eutherian_models": [
            {
                key: row[key] for key in (
                    "species", "protein_accession", "product_class",
                    "transcript_accession", "transcript_evidence_tier",
                    "five_prime_UTR_length_nt",
                )
            }
            for row in non_eutherian
        ],
        "experimental_TSS_support_established": False,
    }
    (RESULTS / "isoform_transcript_support_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
