#!/usr/bin/env python3
"""Audit annotated DNMT3A transcript starts in targeted mammalian genomes.

The audit distinguishes genomic transcript architecture from experimental TSS
evidence. RefSeq mRNA features establish an annotated first exon and TSS, but
predicted XM_/XP_ records do not demonstrate promoter activity.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/promoter_tss"
MANIFEST = RAW / "window_manifest.tsv"
DATA_REPORT = ROOT / "data/raw/ncbi_dataset/data/data_report.jsonl"
PRODUCTS = ROOT / "results/06_isoform_evolution/isoform_product_classification.tsv"
SUPPORT = ROOT / "results/06_isoform_evolution/isoform_transcript_model_support.tsv"
RESULTS = ROOT / "results/06_isoform_evolution"
TARGET_CLASSES = {
    "DNMT3A2_like_downstream_start",
    "leaderless_downstream_core_product",
    "extended_leader_downstream_core_product",
}
PRIMARY_INTERNAL_STARTS = {
    "NM_153759.3",   # human canonical DNMT3A2
    "NM_153743.5",   # mouse annotated isoform 2
    "XM_072633891.1",  # tammar downstream-core model
    "XM_074300067.1",  # dunnart extended-leader model
    "XM_074300068.1",  # dunnart leaderless model
}


def feature_parts(feature) -> list[tuple[int, int]]:
    """Return zero-based, half-open feature intervals."""
    return [(int(part.start), int(part.end)) for part in feature.location.parts]


def tss_local_0based(feature, strand: int) -> int:
    """Return the transcriptional first base in record-local coordinates."""
    parts = feature_parts(feature)
    return min(start for start, _ in parts) if strand == 1 else max(end for _, end in parts) - 1


def first_exon(feature, strand: int) -> tuple[int, int]:
    """Return the first transcribed exon as a zero-based, half-open interval."""
    parts = feature_parts(feature)
    return min(parts) if strand == 1 else max(parts, key=lambda part: part[1])


def genomic_position(fetch_start_1based: int, local_0based: int) -> int:
    return fetch_start_1based + local_0based


def downstream_distance(target_tss: int, outermost_tss: int, strand: int) -> int:
    """Positive values lie downstream in the direction of transcription."""
    return (target_tss - outermost_tss) * strand


def interval_is_exonic(interval: tuple[int, int], exons: list[tuple[int, int]]) -> bool:
    start, end = interval
    return any(exon_start <= start and end <= exon_end for exon_start, exon_end in exons)


def cds_fits_mrna(cds, mrna) -> bool:
    if cds.location.strand != mrna.location.strand:
        return False
    exons = feature_parts(mrna)
    return all(interval_is_exonic(interval, exons) for interval in feature_parts(cds))


def promoter_sequence(record, tss: int, strand: int, upstream: int = 2000, downstream: int = 500):
    """Extract an upstream/downstream sequence in transcriptional orientation."""
    if strand == 1:
        start, end = max(0, tss - upstream), min(len(record), tss + downstream)
        sequence = record.seq[start:end]
    else:
        start, end = max(0, tss - downstream + 1), min(len(record), tss + upstream + 1)
        sequence = record.seq[start:end].reverse_complement()
    return start, end, sequence


def accession_tier(accession: str) -> str:
    if accession.startswith("NM_"):
        return "curated_RefSeq_transcript"
    if accession.startswith("XM_"):
        return "predicted_RefSeq_transcript"
    return "other_transcript"


def load_targets():
    with PRODUCTS.open(newline="") as handle:
        product_rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = {
        row["species"]: row["protein_accession"]
        for row in product_rows
        if row["product_class"] == "selected_full_length"
    }
    product_class = {
        (row["species"], row["protein_accession"]): row["product_class"]
        for row in product_rows
    }
    transcript_products = defaultdict(list)
    with SUPPORT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            transcript_products[(row["species"], row["transcript_accession"])].append(
                (row["protein_accession"], row["product_class"])
            )
    return selected, product_class, transcript_products


def validate_manifest(manifest: list[dict[str, str]]) -> None:
    """Cross-check hand-recorded windows against the packaged NCBI report."""
    reports = {}
    with DATA_REPORT.open() as handle:
        for line in handle:
            report = json.loads(line)
            reports[report["taxname"]] = report
    for item in manifest:
        report = reports[item["species"]]
        if report["geneId"] != item["gene_id"]:
            raise ValueError(f"GeneID mismatch for {item['species']}")
        assembly = next(
            annotation for annotation in report["annotations"]
            if annotation["assemblyAccession"] == item["assembly_accession"]
        )
        if assembly["assemblyName"] != item["assembly"]:
            raise ValueError(f"Assembly-name mismatch for {item['species']}")
        location = next(
            location for location in assembly["genomicLocations"]
            if location["genomicAccessionVersion"] == item["genomic_accession"]
        )
        gene_range = location["genomicRange"]
        expected_strand = 1 if gene_range["orientation"] == "plus" else -1
        if expected_strand != int(item["gene_strand"]):
            raise ValueError(f"Gene-strand mismatch for {item['species']}")
        expected_start = int(gene_range["begin"]) + 1 - 20_000
        expected_stop = int(gene_range["end"]) + 20_000
        if (
            expected_start != int(item["fetch_start_1based"])
            or expected_stop != int(item["fetch_stop_1based"])
        ):
            raise ValueError(f"20-kb window mismatch for {item['species']}")


def main() -> None:
    selected, product_class, transcript_products = load_targets()
    with MANIFEST.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    validate_manifest(manifest)

    rows = []
    promoter_records = []
    primary_promoter_records = []
    primary_bed_rows = []
    checksum_rows = []
    species_summaries = []

    for item in manifest:
        species = item["species"]
        strand = int(item["gene_strand"])
        fetch_start = int(item["fetch_start_1based"])
        path = RAW / item["genbank_file"]
        record = SeqIO.read(path, "genbank")
        checksum_rows.append({
            "file": item["genbank_file"],
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        })

        mrnas = [
            feature for feature in record.features
            if feature.type == "mRNA"
            and feature.qualifiers.get("gene", [""])[0].lower() == "dnmt3a"
        ]
        cdss = [
            feature for feature in record.features
            if feature.type == "CDS"
            and feature.qualifiers.get("gene", [""])[0].lower() == "dnmt3a"
        ]
        if not mrnas:
            raise ValueError(f"No DNMT3A mRNA features found for {species}")
        if any(feature.location.strand != strand for feature in mrnas):
            raise ValueError(f"Manifest/feature strand mismatch for {species}")

        local_tsses = [tss_local_0based(feature, strand) for feature in mrnas]
        outermost_local = min(local_tsses) if strand == 1 else max(local_tsses)
        outermost_genomic = genomic_position(fetch_start, outermost_local)
        selected_protein = selected.get(species, "")
        downstream_target_count = 0

        for mrna in mrnas:
            transcript = mrna.qualifiers.get("transcript_id", [""])[0]
            linked_cdss = [cds for cds in cdss if cds_fits_mrna(cds, mrna)]
            linked_proteins = sorted({
                cds.qualifiers.get("protein_id", [""])[0]
                for cds in linked_cdss
                if cds.qualifiers.get("protein_id", [""])[0]
            })
            supported_products = transcript_products.get((species, transcript), [])
            downstream_products = sorted({
                protein for protein, cls in supported_products if cls in TARGET_CLASSES
            })
            downstream_classes = sorted({
                cls for _, cls in supported_products if cls in TARGET_CLASSES
            })
            is_selected_full = selected_protein in linked_proteins
            is_downstream_target = bool(downstream_products)
            is_primary_target = transcript in PRIMARY_INTERNAL_STARTS
            if is_downstream_target:
                downstream_target_count += 1

            local_tss = tss_local_0based(mrna, strand)
            genomic_tss = genomic_position(fetch_start, local_tss)
            distance = downstream_distance(genomic_tss, outermost_genomic, strand)
            exon_start, exon_end = first_exon(mrna, strand)
            promoter_start, promoter_end, sequence = promoter_sequence(
                record, local_tss, strand
            )
            promoter_id = f"{species.replace(' ', '_')}|{transcript}|TSS_{genomic_tss}"
            promoter_records.append((promoter_id, sequence))
            if is_primary_target:
                primary_promoter_records.append((promoter_id, sequence))
                primary_bed_rows.append((
                    item["genomic_accession"],
                    fetch_start - 1 + promoter_start,
                    fetch_start - 1 + promoter_end,
                    promoter_id,
                    0,
                    "+" if strand == 1 else "-",
                ))

            if distance == 0:
                architecture = "outermost_annotated_start"
            elif distance >= 5000:
                architecture = "distinct_internal_annotated_start"
            else:
                architecture = "proximal_alternative_annotated_start"

            rows.append({
                "species": species,
                "mammal_group": item["mammal_group"],
                "assembly": item["assembly"],
                "assembly_accession": item["assembly_accession"],
                "genomic_accession": item["genomic_accession"],
                "gene_strand": strand,
                "transcript_accession": transcript,
                "transcript_evidence_tier": accession_tier(transcript),
                "annotated_product": mrna.qualifiers.get("product", [""])[0],
                "tss_genomic_1based": genomic_tss,
                "outermost_dnmt3a_tss_genomic_1based": outermost_genomic,
                "distance_downstream_of_outermost_tss_bp": distance,
                "first_exon_genomic_start_1based": genomic_position(fetch_start, exon_start),
                "first_exon_genomic_end_1based": genomic_position(fetch_start, exon_end - 1),
                "first_exon_length_nt": exon_end - exon_start,
                "exon_count": len(feature_parts(mrna)),
                "linked_genbank_protein_accessions": ",".join(linked_proteins),
                "selected_full_length_protein": selected_protein,
                "contains_selected_full_length_CDS": str(is_selected_full).lower(),
                "downstream_core_product_accessions": ",".join(downstream_products),
                "downstream_core_product_classes": ",".join(downstream_classes),
                "is_downstream_core_target": str(is_downstream_target).lower(),
                "is_primary_internal_start_candidate": str(is_primary_target).lower(),
                "start_architecture": architecture,
                "promoter_window_local_start_0based": promoter_start,
                "promoter_window_local_end_0based_exclusive": promoter_end,
                "promoter_window_length_nt": len(sequence),
                "tss_evidence_interpretation": (
                    "curated RefSeq transcript annotation; not a direct TSS assay"
                    if transcript.startswith("NM_")
                    else "predicted RefSeq transcript model; not experimental TSS evidence"
                ),
            })

        species_summaries.append({
            "species": species,
            "annotated_dnmt3a_transcripts": len(mrnas),
            "curated_transcripts": sum(
                feature.qualifiers.get("transcript_id", [""])[0].startswith("NM_")
                for feature in mrnas
            ),
            "downstream_core_target_transcripts": downstream_target_count,
            "has_distinct_internal_downstream_core_model": any(
                row["species"] == species
                and row["is_downstream_core_target"] == "true"
                and row["start_architecture"] == "distinct_internal_annotated_start"
                for row in rows
            ),
        })

    fields = list(rows[0])
    with (RESULTS / "target_transcript_tss.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    with (RESULTS / "target_promoter_windows.fasta").open("w") as handle:
        for identifier, sequence in promoter_records:
            handle.write(f">{identifier}\n")
            text = str(sequence).upper()
            for start in range(0, len(text), 80):
                handle.write(text[start:start + 80] + "\n")

    with (RESULTS / "primary_internal_start_promoters.fasta").open("w") as handle:
        for identifier, sequence in primary_promoter_records:
            handle.write(f">{identifier}\n")
            text = str(sequence).upper()
            for start in range(0, len(text), 80):
                handle.write(text[start:start + 80] + "\n")

    with (RESULTS / "primary_internal_start_promoters.bed").open("w") as handle:
        for row in primary_bed_rows:
            handle.write("\t".join(map(str, row)) + "\n")

    with (RAW / "SHA256SUMS.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(checksum_rows[0]))
        writer.writeheader()
        writer.writerows(checksum_rows)

    target_rows = [row for row in rows if row["is_downstream_core_target"] == "true"]
    summary = {
        "species_audited": len(manifest),
        "annotated_dnmt3a_transcripts": len(rows),
        "downstream_core_target_transcripts": len(target_rows),
        "distinct_internal_downstream_core_transcripts": sum(
            row["start_architecture"] == "distinct_internal_annotated_start"
            for row in target_rows
        ),
        "species": species_summaries,
        "target_transcripts": [
            {
                key: row[key] for key in (
                    "species", "transcript_accession", "transcript_evidence_tier",
                    "tss_genomic_1based", "distance_downstream_of_outermost_tss_bp",
                    "first_exon_length_nt", "downstream_core_product_accessions",
                    "downstream_core_product_classes", "start_architecture",
                    "tss_evidence_interpretation",
                )
            }
            for row in target_rows
        ],
        "primary_internal_start_transcripts": sorted(PRIMARY_INTERNAL_STARTS),
        "primary_internal_start_export_count": len(primary_promoter_records),
        "experimental_TSS_support_established": False,
        "inference_gate": (
            "Annotated internal starts support promoter-architecture candidates; "
            "evolutionary gain/loss requires orthologous sequence and direct "
            "transcription-start evidence."
        ),
    }
    (RESULTS / "target_tss_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
