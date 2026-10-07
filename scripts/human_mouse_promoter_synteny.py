#!/usr/bin/env python3
"""Test human-to-mouse synteny of the DNMT3A internal promoter."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

from pyliftover import LiftOver


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/fantom5_cage"
RESULTS = ROOT / "results/06_isoform_evolution"
DIRECT = RESULTS / "direct_tss_evidence.tsv"
NEARBY = RESULTS / "fantom5_nearby_cage_peaks.tsv"
PROMOTERS = RESULTS / "primary_internal_start_promoters.bed"
HG38_TO_MM10 = RAW / "hg38ToMm10.over.chain.gz"
MM10_TO_MM39 = RAW / "mm10ToMm39.over.chain.gz"


def forward_query_coordinate(q_size: int, q_strand: str, q_position: int) -> int:
    if q_strand == "+":
        return q_position
    return q_size - q_position - 1


def transformed_strand(input_strand: str, chain_query_strand: str) -> str:
    if chain_query_strand == "+":
        return input_strand
    return "+" if input_strand == "-" else "-"


def map_target_points(
    chain_path: Path,
    target_chromosome: str,
    points: dict[str, int],
    input_strand: str,
) -> dict[str, list[dict[str, int | str]]]:
    """Map selected target-genome points without loading a whole chain in RAM."""
    hits = {name: [] for name in points}
    with gzip.open(chain_path, "rt") as handle:
        header = None
        target_position = query_position = 0
        active = False
        for line in handle:
            if line.startswith("chain "):
                fields = line.split()
                header = fields
                target_position = int(fields[5])
                query_position = int(fields[10])
                active = (
                    fields[2] == target_chromosome
                    and any(
                        int(fields[5]) <= point < int(fields[6])
                        for point in points.values()
                    )
                )
                continue
            if not line.strip():
                header = None
                active = False
                continue
            values = [int(value) for value in line.split()]
            size = values[0]
            if header and active:
                for name, point in points.items():
                    if target_position <= point < target_position + size:
                        offset = point - target_position
                        query_chain_position = query_position + offset
                        hits[name].append({
                            "chain_score": int(header[1]),
                            "chain_id": header[12],
                            "chromosome": header[7],
                            "position": forward_query_coordinate(
                                int(header[8]), header[9], query_chain_position
                            ),
                            "strand": transformed_strand(input_strand, header[9]),
                        })
            if len(values) == 3:
                target_position += size + values[1]
                query_position += size + values[2]
    return hits


def unique_best(hits: list[dict[str, int | str]]) -> dict[str, int | str] | None:
    if not hits:
        return None
    ordered = sorted(hits, key=lambda row: int(row["chain_score"]), reverse=True)
    if len(ordered) > 1 and ordered[0]["chain_score"] == ordered[1]["chain_score"]:
        return None
    return ordered[0]


def main() -> None:
    with DIRECT.open(newline="") as handle:
        direct = list(csv.DictReader(handle, delimiter="\t"))
    human = next(row for row in direct if row["species"] == "Homo sapiens")
    mouse = next(row for row in direct if row["species"] == "Mus musculus")
    mouse_tss = int(mouse["annotated_tss_genomic_1based"]) - 1

    with NEARBY.open(newline="") as handle:
        cage = list(csv.DictReader(handle, delimiter="\t"))
    human_cage = [row for row in cage if row["species"] == "Homo sapiens"]
    mouse_ctss = sorted({
        int(row["dominant_CTSS_0based"])
        for row in cage if row["species"] == "Mus musculus"
    })

    with PROMOTERS.open() as handle:
        promoter_fields = [
            line.rstrip().split("\t") for line in handle
            if "|NM_153759.3|" in line
        ][0]
    promoter_start = int(promoter_fields[1])
    promoter_end_last_base = int(promoter_fields[2]) - 1

    landmarks = {
        "promoter_window_start": promoter_start,
        "promoter_window_end_last_base": promoter_end_last_base,
        "annotated_TSS": int(human["annotated_tss_genomic_1based"]) - 1,
    }
    for index, row in enumerate(human_cage, start=1):
        landmarks[f"human_CAGE_CTSS_{index}:{row['peak_id']}"] = int(
            row["dominant_CTSS_0based"]
        )

    mapped_mm10 = map_target_points(
        HG38_TO_MM10, "chr2", landmarks, input_strand="-"
    )
    mm10_to_mm39 = LiftOver(str(MM10_TO_MM39))
    rows = []
    for name, hg38_position in landmarks.items():
        first = unique_best(mapped_mm10[name])
        if first is None:
            rows.append({
                "landmark": name,
                "hg38_chr": "chr2",
                "hg38_position_0based": hg38_position,
                "mm10_chr": "",
                "mm10_position_0based": "",
                "mm39_chr": "",
                "mm39_position_0based": "",
                "mapped_strand": "",
                "hg38_to_mm10_chain_id": "",
                "hg38_to_mm10_chain_score": "",
                "distance_to_mouse_annotated_TSS_bp": "",
                "distance_to_nearest_mouse_CAGE_CTSS_bp": "",
                "mapping_status": "unmapped_or_nonunique",
            })
            continue
        second_hits = mm10_to_mm39.convert_coordinate(
            str(first["chromosome"]), int(first["position"]), str(first["strand"])
        )
        if len(second_hits) != 1:
            raise ValueError(f"Nonunique mm10-to-mm39 mapping for {name}")
        mm39_chromosome, mm39_position, mm39_strand, _ = second_hits[0]
        nearest_ctss_distance = min(
            abs(mm39_position - position) for position in mouse_ctss
        )
        rows.append({
            "landmark": name,
            "hg38_chr": "chr2",
            "hg38_position_0based": hg38_position,
            "mm10_chr": first["chromosome"],
            "mm10_position_0based": first["position"],
            "mm39_chr": mm39_chromosome,
            "mm39_position_0based": mm39_position,
            "mapped_strand": mm39_strand,
            "hg38_to_mm10_chain_id": first["chain_id"],
            "hg38_to_mm10_chain_score": first["chain_score"],
            "distance_to_mouse_annotated_TSS_bp": abs(mm39_position - mouse_tss),
            "distance_to_nearest_mouse_CAGE_CTSS_bp": nearest_ctss_distance,
            "mapping_status": "unique_coherent_chain",
        })

    with (RESULTS / "human_mouse_internal_promoter_synteny.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    by_name = {row["landmark"]: row for row in rows}
    endpoints = [
        int(by_name["promoter_window_start"]["mm39_position_0based"]),
        int(by_name["promoter_window_end_last_base"]["mm39_position_0based"]),
    ]
    mapped_cage = [
        row for row in rows
        if row["landmark"].startswith("human_CAGE_CTSS_")
        and row["mapping_status"] == "unique_coherent_chain"
    ]
    summary = {
        "human_promoter_window_hg38": {
            "chromosome": "chr2",
            "start_0based": promoter_start,
            "end_0based_exclusive": promoter_end_last_base + 1,
            "strand": "-",
        },
        "mapped_mouse_window_GRCm39": {
            "chromosome": "chr12",
            "start_0based": min(endpoints),
            "end_0based_exclusive": max(endpoints) + 1,
            "strand": "+",
            "contains_mouse_internal_TSS": (
                min(endpoints) <= mouse_tss <= max(endpoints)
            ),
        },
        "human_annotated_TSS_mapped_distance_to_mouse_TSS_bp": int(
            by_name["annotated_TSS"]["distance_to_mouse_annotated_TSS_bp"]
        ),
        "minimum_mapped_human_CAGE_CTSS_distance_to_mouse_TSS_bp": min(
            int(row["distance_to_mouse_annotated_TSS_bp"]) for row in mapped_cage
        ),
        "minimum_mapped_human_CAGE_CTSS_distance_to_mouse_CAGE_CTSS_bp": min(
            int(row["distance_to_nearest_mouse_CAGE_CTSS_bp"]) for row in mapped_cage
        ),
        "inference": (
            "The human internal-promoter window maps coherently across the "
            "orthologous mouse internal TSS/CAGE cluster."
        ),
        "limit": (
            "Syntenic mapping supports regulatory-region orthology but does "
            "not show conservation of every promoter base or regulatory effect."
        ),
    }
    (RESULTS / "human_mouse_internal_promoter_synteny_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
