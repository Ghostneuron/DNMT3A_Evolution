#!/usr/bin/env python3
"""Compare DNMT3A internal starts relative to conserved coding-exon anchors."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from Bio import SeqIO

from promoter_tss_audit import feature_parts, first_exon, tss_local_0based


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/promoter_tss"
MANIFEST = RAW / "window_manifest.tsv"
PRODUCTS = ROOT / "results/06_isoform_evolution/isoform_product_classification.tsv"
SUPPORT = ROOT / "results/06_isoform_evolution/isoform_transcript_model_support.tsv"
TARGETS = ROOT / "results/06_isoform_evolution/target_transcript_tss.tsv"
RESULTS = ROOT / "results/06_isoform_evolution"
DEFAULT_CORE_AA = 212


def ordered_parts(feature, strand: int) -> list[tuple[int, int]]:
    return sorted(
        feature_parts(feature),
        key=lambda interval: interval[0],
        reverse=strand == -1,
    )


def overlaps(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return max(left[0], right[0]) < min(left[1], right[1])


def transcriptional_start(interval: tuple[int, int], strand: int) -> int:
    return interval[0] if strand == 1 else interval[1] - 1


def spliced_bases_before_coordinate(feature, coordinate: int, strand: int) -> tuple[int, int]:
    total = 0
    for index, (start, end) in enumerate(ordered_parts(feature, strand)):
        if start <= coordinate < end:
            total += coordinate - start if strand == 1 else end - 1 - coordinate
            return total, index
        total += end - start
    raise ValueError(f"Coordinate {coordinate} is not in feature")


def coordinate_at_coding_offset(cds, nt_offset: int, strand: int) -> int:
    remaining = nt_offset
    for start, end in ordered_parts(cds, strand):
        length = end - start
        if remaining < length:
            return start + remaining if strand == 1 else end - 1 - remaining
        remaining -= length
    raise ValueError("Coding offset exceeds CDS")


def main() -> None:
    with MANIFEST.open(newline="") as handle:
        manifest = {row["species"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with PRODUCTS.open(newline="") as handle:
        products = list(csv.DictReader(handle, delimiter="\t"))
    selected = {
        row["species"]: row["protein_accession"]
        for row in products if row["product_class"] == "selected_full_length"
    }
    product_rows = {
        (row["species"], row["protein_accession"]): row for row in products
    }
    with SUPPORT.open(newline="") as handle:
        support = list(csv.DictReader(handle, delimiter="\t"))
    supported_products = {}
    for row in support:
        supported_products.setdefault(
            (row["species"], row["transcript_accession"]), []
        ).append(row["protein_accession"])
    with TARGETS.open(newline="") as handle:
        primary = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["is_primary_internal_start_candidate"] == "true"
        ]

    rows = []
    for target in primary:
        species = target["species"]
        item = manifest[species]
        strand = int(item["gene_strand"])
        record = SeqIO.read(RAW / item["genbank_file"], "genbank")
        mrnas = {
            feature.qualifiers.get("transcript_id", [""])[0]: feature
            for feature in record.features if feature.type == "mRNA"
        }
        cdss = {
            feature.qualifiers.get("protein_id", [""])[0]: feature
            for feature in record.features if feature.type == "CDS"
        }
        transcript = target["transcript_accession"]
        target_protein = next(
            protein for protein in supported_products[(species, transcript)]
            if protein in cdss
        )
        full_protein = selected[species]
        full_cds = cdss[full_protein]
        target_cds = cdss[target_protein]
        target_parts = ordered_parts(target_cds, strand)
        full_parts = ordered_parts(full_cds, strand)
        shared_index = next(
            index for index, part in enumerate(target_parts)
            if any(overlaps(part, full_part) for full_part in full_parts)
        )
        shared_part = target_parts[shared_index]
        anchor = transcriptional_start(shared_part, strand)
        mrna = mrnas[transcript]
        tss = tss_local_0based(mrna, strand)
        spliced_before, exon_index = spliced_bases_before_coordinate(
            mrna, anchor, strand
        )
        full_nt_before, _ = spliced_bases_before_coordinate(
            full_cds, anchor, strand
        )
        target_nt_before, _ = spliced_bases_before_coordinate(
            target_cds, anchor, strand
        )
        classification = product_rows[(species, target_protein)]
        fetch_start = int(item["fetch_start_1based"])
        rows.append({
            "species": species,
            "mammal_group": item["mammal_group"],
            "transcript_accession": transcript,
            "target_protein_accession": target_protein,
            "target_product_class": classification["product_class"],
            "selected_full_protein_accession": full_protein,
            "gene_strand": strand,
            "annotated_tss_genomic_1based": fetch_start + tss,
            "shared_core_anchor_genomic_1based": fetch_start + anchor,
            "genomic_distance_tss_to_shared_anchor_bp": (anchor - tss) * strand,
            "spliced_transcript_bases_before_shared_anchor": spliced_before,
            "transcript_exons_before_shared_anchor": exon_index,
            "target_CDS_segments_before_shared_anchor": shared_index,
            "target_coding_nt_before_shared_anchor": target_nt_before,
            "full_CDS_nt_before_shared_anchor": full_nt_before,
            "full_protein_anchor_aa_0based": full_nt_before // 3,
            "sequence_classification_core_start_aa_0based": classification[
                "shared_block_full_start_0based"
            ],
            "architecture_interpretation": (
                "Internal transcript start followed by a coding segment shared "
                "with the full-length DNMT3A product"
            ),
            "evidence_tier": (
                "curated_transcript_architecture"
                if transcript.startswith("NM_")
                else "predicted_transcript_architecture"
            ),
        })

    # Platypus has no primary internal-start model. Locate the homologous
    # full-length coding neighborhood for an explicit annotation-absence row.
    species = "Ornithorhynchus anatinus"
    item = manifest[species]
    strand = int(item["gene_strand"])
    record = SeqIO.read(RAW / item["genbank_file"], "genbank")
    full_protein = selected[species]
    full_cds = next(
        feature for feature in record.features
        if feature.type == "CDS"
        and feature.qualifiers.get("protein_id", [""])[0] == full_protein
    )
    anchor = coordinate_at_coding_offset(full_cds, DEFAULT_CORE_AA * 3, strand)
    fetch_start = int(item["fetch_start_1based"])
    rows.append({
        "species": species,
        "mammal_group": item["mammal_group"],
        "transcript_accession": "",
        "target_protein_accession": "",
        "target_product_class": "no_internal_start_model_in_current_annotation",
        "selected_full_protein_accession": full_protein,
        "gene_strand": strand,
        "annotated_tss_genomic_1based": "",
        "shared_core_anchor_genomic_1based": fetch_start + anchor,
        "genomic_distance_tss_to_shared_anchor_bp": "",
        "spliced_transcript_bases_before_shared_anchor": "",
        "transcript_exons_before_shared_anchor": "",
        "target_CDS_segments_before_shared_anchor": "",
        "target_coding_nt_before_shared_anchor": "",
        "full_CDS_nt_before_shared_anchor": DEFAULT_CORE_AA * 3,
        "full_protein_anchor_aa_0based": DEFAULT_CORE_AA,
        "sequence_classification_core_start_aa_0based": "",
        "architecture_interpretation": (
            "Homologous full-length coding neighborhood present; no distinct "
            "internal transcript start is annotated"
        ),
        "evidence_tier": "annotation_absence_not_biological_absence",
    })

    with (RESULTS / "internal_start_exon_architecture.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    target_rows = [row for row in rows if row["transcript_accession"]]
    summary = {
        "primary_internal_start_models": len(target_rows),
        "species_with_internal_start_architecture": sorted({
            row["species"] for row in target_rows
        }),
        "curated_models": sum(
            row["evidence_tier"] == "curated_transcript_architecture"
            for row in target_rows
        ),
        "predicted_models": sum(
            row["evidence_tier"] == "predicted_transcript_architecture"
            for row in target_rows
        ),
        "full_protein_anchor_aa_range_0based": [
            min(int(row["full_protein_anchor_aa_0based"]) for row in target_rows),
            max(int(row["full_protein_anchor_aa_0based"]) for row in target_rows),
        ],
        "platypus_annotation_status": (
            "Homologous core coding neighborhood present; no internal-start "
            "transcript in current RefSeq annotation"
        ),
        "inference": (
            "Human, mouse, tammar, and dunnart internal transcripts converge "
            "on the same conserved DNMT3A coding neighborhood."
        ),
        "limit": (
            "The marsupial models are predicted; shared architecture alone "
            "does not demonstrate active promoters."
        ),
    }
    (RESULTS / "internal_start_exon_architecture_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
