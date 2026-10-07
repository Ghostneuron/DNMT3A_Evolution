#!/usr/bin/env python3
"""Overlap primary DNMT3A internal starts with direct FANTOM5 CAGE peaks."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

from pyliftover import LiftOver


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/fantom5_cage"
MANIFEST = RAW / "source_manifest.tsv"
TARGETS = ROOT / "results/06_isoform_evolution/target_transcript_tss.tsv"
RESULTS = ROOT / "results/06_isoform_evolution"
SEARCH_RADIUS = 500
CHROMOSOMES = {
    ("Homo sapiens", "NC_000002.12"): "chr2",
    ("Mus musculus", "NC_000078.7"): "chr12",
}


def point_interval_distance(point: int, start: int, end: int) -> int:
    """Distance from a zero-based point to a half-open BED interval."""
    if point < start:
        return start - point
    if point >= end:
        return point - end + 1
    return 0


def evidence_tier(interval_distance: int, dominant_distance: int) -> str:
    if interval_distance == 0 and dominant_distance <= 5:
        return "direct_CAGE_peak_exact_or_near_exact"
    if interval_distance <= 10 or dominant_distance <= 50:
        return "direct_CAGE_peak_near_annotated_TSS"
    if interval_distance <= SEARCH_RADIUS:
        return "direct_CAGE_peak_within_500bp"
    return "no_CAGE_peak_within_500bp"


def validate_files() -> dict[str, dict[str, str]]:
    with MANIFEST.open(newline="") as handle:
        sources = {
            row["source_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    for row in sources.values():
        path = RAW / row["local_file"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != row["sha256"]:
            raise ValueError(f"SHA-256 mismatch for {path}")
    return sources


def parse_bed(path: Path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 8:
                continue
            yield {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "peak_id": fields[3],
                "score": fields[4],
                "strand": fields[5],
                "dominant_ctss": int(fields[6]),
            }


def lift_peak(peak: dict, lifter: LiftOver) -> dict | None:
    """Lift a BED interval and dominant CTSS, requiring a unique coherent map."""
    mapped = []
    for position in (peak["start"], peak["end"] - 1, peak["dominant_ctss"]):
        hits = lifter.convert_coordinate(peak["chrom"], position, peak["strand"])
        if len(hits) != 1:
            return None
        mapped.append(hits[0])
    chromosomes = {hit[0] for hit in mapped}
    strands = {hit[2] for hit in mapped}
    if len(chromosomes) != 1 or len(strands) != 1:
        return None
    endpoints = [mapped[0][1], mapped[1][1]]
    lifted = dict(peak)
    lifted.update({
        "chrom": mapped[0][0],
        "start": min(endpoints),
        "end": max(endpoints) + 1,
        "strand": mapped[0][2],
        "dominant_ctss": mapped[2][1],
    })
    return lifted


def load_primary_targets() -> list[dict[str, str]]:
    with TARGETS.open(newline="") as handle:
        return [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["is_primary_internal_start_candidate"] == "true"
        ]


def main() -> None:
    sources = validate_files()
    targets = load_primary_targets()
    assay_targets = [
        row for row in targets if row["species"] in {"Homo sapiens", "Mus musculus"}
    ]
    target_index = {
        (row["species"], CHROMOSOMES[(row["species"], row["genomic_accession"])]): row
        for row in assay_targets
    }

    nearby = []
    human_path = RAW / sources["fantom5_hg38_fair_new_peaks"]["local_file"]
    for peak in parse_bed(human_path):
        target = target_index.get(("Homo sapiens", peak["chrom"]))
        if target:
            collect_peak(target, peak, "GRCh38/hg38", False, nearby)

    mouse_path = RAW / sources["fantom5_mm10_fair_new_peaks"]["local_file"]
    chain_path = RAW / sources["ucsc_mm10_to_mm39_chain"]["local_file"]
    lifter = LiftOver(str(chain_path))
    for source_peak in parse_bed(mouse_path):
        if source_peak["chrom"] != "chr12":
            continue
        peak = lift_peak(source_peak, lifter)
        target = target_index.get(("Mus musculus", peak["chrom"])) if peak else None
        if target:
            collect_peak(target, peak, "GRCm39", True, nearby)

    nearby.sort(key=lambda row: (
        row["species"], row["transcript_accession"],
        int(row["distance_to_peak_interval_bp"]),
        int(row["distance_to_dominant_CTSS_bp"]),
    ))
    nearby_fields = list(nearby[0])
    with (RESULTS / "fantom5_nearby_cage_peaks.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=nearby_fields)
        writer.writeheader()
        writer.writerows(nearby)

    summaries = []
    for target in targets:
        matches = [
            row for row in nearby
            if row["species"] == target["species"]
            and row["transcript_accession"] == target["transcript_accession"]
        ]
        best = matches[0] if matches else None
        if target["species"] in {"Homo sapiens", "Mus musculus"}:
            status = (
                best["evidence_tier"]
                if best else "no_CAGE_peak_within_500bp"
            )
            source = "FANTOM5 CAGE"
        else:
            status = "not_assessed_no_species_specific_direct_TSS_dataset"
            source = ""
        summaries.append({
            "species": target["species"],
            "mammal_group": target["mammal_group"],
            "transcript_accession": target["transcript_accession"],
            "transcript_evidence_tier": target["transcript_evidence_tier"],
            "annotated_tss_genomic_1based": target["tss_genomic_1based"],
            "analysis_assembly": target["assembly"],
            "direct_tss_source": source,
            "direct_tss_evidence_status": status,
            "closest_peak_id": best["peak_id"] if best else "",
            "closest_peak_interval_distance_bp": (
                best["distance_to_peak_interval_bp"] if best else ""
            ),
            "closest_dominant_CTSS_distance_bp": (
                best["distance_to_dominant_CTSS_bp"] if best else ""
            ),
            "same_strand": best["same_strand"] if best else "",
            "interpretation": (
                "Direct capped-RNA initiation evidence near the annotated internal start"
                if best else
                "Annotation model only; absence was not tested with a compatible direct-TSS atlas"
            ),
        })

    with (RESULTS / "direct_tss_evidence.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(summaries[0])
        )
        writer.writeheader()
        writer.writerows(summaries)

    summary = {
        "primary_internal_start_models": len(targets),
        "models_with_compatible_FANTOM5_assay": len(assay_targets),
        "models_with_direct_CAGE_peak_within_500bp": sum(
            row["direct_tss_evidence_status"].startswith("direct_CAGE")
            for row in summaries
        ),
        "annotation_only_models_without_compatible_direct_TSS_atlas": sum(
            row["direct_tss_evidence_status"].startswith("not_assessed")
            for row in summaries
        ),
        "results": summaries,
        "interpretation_limit": (
            "FANTOM5 establishes transcription initiation in human and mouse "
            "near these loci; it does not validate marsupial models or by itself "
            "establish promoter orthology."
        ),
    }
    (RESULTS / "direct_tss_evidence_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


def collect_peak(
    target: dict[str, str],
    peak: dict,
    analysis_assembly: str,
    lifted: bool,
    output: list[dict[str, str]],
) -> None:
    target_tss = int(target["tss_genomic_1based"]) - 1
    interval_distance = point_interval_distance(
        target_tss, peak["start"], peak["end"]
    )
    if interval_distance > SEARCH_RADIUS:
        return
    dominant_distance = abs(target_tss - peak["dominant_ctss"])
    target_strand = "+" if target["gene_strand"] == "1" else "-"
    same_strand = peak["strand"] == target_strand
    if not same_strand:
        return
    output.append({
        "species": target["species"],
        "transcript_accession": target["transcript_accession"],
        "analysis_assembly": analysis_assembly,
        "chromosome": peak["chrom"],
        "annotated_tss_0based": target_tss,
        "peak_start_0based": peak["start"],
        "peak_end_0based_exclusive": peak["end"],
        "dominant_CTSS_0based": peak["dominant_ctss"],
        "peak_id": peak["peak_id"],
        "peak_score": peak["score"],
        "peak_strand": peak["strand"],
        "same_strand": str(same_strand).lower(),
        "distance_to_peak_interval_bp": interval_distance,
        "distance_to_dominant_CTSS_bp": dominant_distance,
        "mouse_mm10_to_mm39_lifted": str(lifted).lower(),
        "evidence_tier": evidence_tier(interval_distance, dominant_distance),
    })


if __name__ == "__main__":
    main()
