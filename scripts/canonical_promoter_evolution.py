#!/usr/bin/env python3
"""Analyze canonical full-length DNMT3A1 promoter regions across mammals."""

from __future__ import annotations

import csv
import json
from difflib import SequenceMatcher
from itertools import combinations
from pathlib import Path

from Bio import SeqIO
from pyliftover import LiftOver

try:
    from scripts.direct_tss_evidence import (
        collect_peak,
        lift_peak,
        parse_bed,
        validate_files,
    )
    from scripts.human_mouse_promoter_synteny import (
        map_target_points,
        unique_best,
    )
    from scripts.promoter_sequence_similarity import kmers
except ModuleNotFoundError:
    from direct_tss_evidence import (
        collect_peak,
        lift_peak,
        parse_bed,
        validate_files,
    )
    from human_mouse_promoter_synteny import (
        map_target_points,
        unique_best,
    )
    from promoter_sequence_similarity import kmers


ROOT = Path(__file__).resolve().parents[1]
RAW_CAGE = ROOT / "data/raw/fantom5_cage"
RESULTS = ROOT / "results/06_isoform_evolution"
TARGETS = RESULTS / "target_transcript_tss.tsv"
ALL_PROMOTERS = RESULTS / "target_promoter_windows.fasta"
OUT = RESULTS / "canonical_dnmt3a1_promoter"
HG38_TO_MM10 = RAW_CAGE / "hg38ToMm10.over.chain.gz"
MM10_TO_MM39 = RAW_CAGE / "mm10ToMm39.over.chain.gz"
PRIMARY_FULL_LENGTH = {
    "Homo sapiens": "NM_175629.2",
    "Mus musculus": "NM_007872.5",
    "Notamacropus eugenii": "XM_072633885.1",
    "Ornithorhynchus anatinus": "XM_029071763.1",
    "Sminthopsis crassicaudata": "XM_074300062.1",
}
CURATED_FULL_LENGTH_TSS_CLUSTER = {
    "Homo sapiens": {"NM_175629.2", "NM_022552.5"},
    "Mus musculus": {
        "NM_007872.5",
        "NM_001271753.2",
        "NM_001421857.1",
    },
}
CHROMOSOMES = {
    "Homo sapiens": "chr2",
    "Mus musculus": "chr12",
}
KMER_SIZES = (11, 15, 19)
REGIONS = {
    "upstream_2000bp": (0, 2000),
    "downstream_500bp": (2000, 2500),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty output: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def select_primary_targets(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected = []
    for species, accession in PRIMARY_FULL_LENGTH.items():
        matches = [
            row for row in rows
            if row["species"] == species
            and row["transcript_accession"] == accession
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one canonical target for {species} {accession}; "
                f"found {len(matches)}"
            )
        row = matches[0]
        if row["contains_selected_full_length_CDS"] != "true":
            raise ValueError(f"{species} canonical target is not full length")
        selected.append(row)
    return selected


def select_tss_cluster(
    rows: list[dict[str, str]],
    primary: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Select curated eutherian DNMT3A1 starts plus non-model representatives."""
    selected = []
    for species, accessions in CURATED_FULL_LENGTH_TSS_CLUSTER.items():
        matches = [
            row for row in rows
            if row["species"] == species
            and row["transcript_accession"] in accessions
        ]
        if {row["transcript_accession"] for row in matches} != accessions:
            raise ValueError(f"Incomplete curated DNMT3A1 TSS cluster: {species}")
        selected.extend(matches)
    selected.extend(
        row for row in primary
        if row["species"] not in CURATED_FULL_LENGTH_TSS_CLUSTER
    )
    if any(row["contains_selected_full_length_CDS"] != "true" for row in selected):
        raise ValueError("TSS cluster contains a non-full-length transcript")
    return selected


def promoter_records(
    targets: list[dict[str, str]],
) -> dict[str, tuple[dict[str, str], str]]:
    all_records = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(ALL_PROMOTERS, "fasta")
    }
    records = {}
    for target in targets:
        prefix = (
            f"{target['species'].replace(' ', '_')}|"
            f"{target['transcript_accession']}|"
        )
        matches = [
            (record_id, sequence)
            for record_id, sequence in all_records.items()
            if record_id.startswith(prefix)
        ]
        if len(matches) != 1:
            raise ValueError(f"Expected one promoter sequence for {prefix}")
        record_id, sequence = matches[0]
        if len(sequence) != 2500:
            raise ValueError(f"Unexpected promoter length for {record_id}")
        records[record_id] = (target, sequence)
    return records


def sequence_similarity(
    records: dict[str, tuple[dict[str, str], str]],
) -> list[dict[str, object]]:
    rows = []
    for (left_id, (_, left)), (right_id, (_, right)) in combinations(
        records.items(), 2
    ):
        for region, (start, end) in REGIONS.items():
            left_region = left[start:end]
            right_region = right[start:end]
            longest = SequenceMatcher(
                None, left_region, right_region, autojunk=False
            ).find_longest_match()
            for size in KMER_SIZES:
                left_kmers = kmers(left_region, size)
                right_kmers = kmers(right_region, size)
                shared = left_kmers & right_kmers
                union = left_kmers | right_kmers
                rows.append({
                    "left_model": left_id,
                    "right_model": right_id,
                    "region": region,
                    "kmer_size": size,
                    "left_unique_kmers": len(left_kmers),
                    "right_unique_kmers": len(right_kmers),
                    "shared_unique_kmers": len(shared),
                    "kmer_jaccard": len(shared) / len(union) if union else 0,
                    "longest_exact_block_length": longest.size,
                    "longest_exact_block_left_offset_in_region": longest.a,
                    "longest_exact_block_right_offset_in_region": longest.b,
                    "longest_exact_block_sequence": left_region[
                        longest.a:longest.a + longest.size
                    ],
                    "interpretation_limit": (
                        "Descriptive locus-anchored exact similarity; it is "
                        "not a promoter-function or orthology test"
                    ),
                })
    return rows


def collect_cage(
    targets: list[dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    sources = validate_files()
    by_species = {
        species: [row for row in targets if row["species"] == species]
        for species in CHROMOSOMES
    }
    nearby: list[dict[str, str]] = []
    human_path = (
        RAW_CAGE / sources["fantom5_hg38_fair_new_peaks"]["local_file"]
    )
    for peak in parse_bed(human_path):
        if peak["chrom"] == "chr2":
            for target in by_species["Homo sapiens"]:
                collect_peak(
                    target, peak, "GRCh38/hg38", False, nearby
                )

    mouse_path = (
        RAW_CAGE / sources["fantom5_mm10_fair_new_peaks"]["local_file"]
    )
    lifter = LiftOver(
        str(RAW_CAGE / sources["ucsc_mm10_to_mm39_chain"]["local_file"])
    )
    for source_peak in parse_bed(mouse_path):
        if source_peak["chrom"] != "chr12":
            continue
        peak = lift_peak(source_peak, lifter)
        if peak:
            for target in by_species["Mus musculus"]:
                collect_peak(
                    target, peak, "GRCm39", True, nearby
                )
    nearby.sort(
        key=lambda row: (
            row["species"],
            int(row["distance_to_peak_interval_bp"]),
            int(row["distance_to_dominant_CTSS_bp"]),
        )
    )

    summaries = []
    for target in targets:
        matches = [
            row for row in nearby
            if row["species"] == target["species"]
            and row["transcript_accession"] == target["transcript_accession"]
        ]
        best = matches[0] if matches else None
        assayed = target["species"] in CHROMOSOMES
        summaries.append({
            "species": target["species"],
            "transcript_accession": target["transcript_accession"],
            "transcript_evidence_tier": target["transcript_evidence_tier"],
            "annotated_tss_genomic_1based": target["tss_genomic_1based"],
            "direct_tss_source": "FANTOM5 CAGE" if assayed else "",
            "direct_tss_evidence_status": (
                best["evidence_tier"]
                if best else (
                    "no_CAGE_peak_within_500bp"
                    if assayed else
                    "not_assessed_no_species_specific_direct_TSS_dataset"
                )
            ),
            "closest_peak_id": best["peak_id"] if best else "",
            "closest_peak_interval_distance_bp": (
                best["distance_to_peak_interval_bp"] if best else ""
            ),
            "closest_dominant_CTSS_distance_bp": (
                best["distance_to_dominant_CTSS_bp"] if best else ""
            ),
            "same_strand": best["same_strand"] if best else "",
            "interpretation": (
                "Direct capped-RNA initiation evidence near the canonical "
                "full-length transcript start"
                if best else (
                    "No compatible direct-TSS dataset was assessed"
                    if not assayed else
                    "No same-strand FANTOM5 peak was found within 500 bp"
                )
            ),
        })
    return nearby, summaries


def human_mouse_synteny(
    targets: list[dict[str, str]],
    nearby: list[dict[str, object]],
    records: dict[str, tuple[dict[str, str], str]],
) -> tuple[list[dict[str, object]], dict]:
    human_targets = [
        row for row in targets if row["species"] == "Homo sapiens"
    ]
    mouse_targets = [
        row for row in targets if row["species"] == "Mus musculus"
    ]
    human = next(
        row for row in human_targets
        if row["transcript_accession"] == PRIMARY_FULL_LENGTH["Homo sapiens"]
    )
    human_tss = int(human["tss_genomic_1based"]) - 1
    mouse_tsses = sorted({
        int(row["tss_genomic_1based"]) - 1 for row in mouse_targets
    })
    human_record = next(
        record_id for record_id in records if record_id.startswith("Homo_sapiens")
    )
    human_target = records[human_record][0]
    local_start = int(human_target["promoter_window_local_start_0based"])
    local_end = int(human_target["promoter_window_local_end_0based_exclusive"])
    promoter_start = human_tss - 499
    promoter_end_last_base = human_tss + 2000
    if local_end - local_start != 2500:
        raise ValueError("Canonical human promoter window is not 2500 bp")

    human_cage = [
        row for row in nearby if row["species"] == "Homo sapiens"
    ]
    mouse_ctss = sorted({
        int(row["dominant_CTSS_0based"])
        for row in nearby if row["species"] == "Mus musculus"
    })
    landmarks = {
        "promoter_window_start": promoter_start,
        "promoter_window_end_last_base": promoter_end_last_base,
    }
    for target in human_targets:
        landmarks[
            f"human_annotated_TSS:{target['transcript_accession']}"
        ] = int(target["tss_genomic_1based"]) - 1
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
                "hg38_position_0based": hg38_position,
                "mm39_position_0based": "",
                "mapped_strand": "",
                "distance_to_nearest_mouse_annotated_TSS_bp": "",
                "distance_to_nearest_mouse_CAGE_CTSS_bp": "",
                "mapping_status": "unmapped_or_nonunique",
            })
            continue
        second_hits = mm10_to_mm39.convert_coordinate(
            str(first["chromosome"]),
            int(first["position"]),
            str(first["strand"]),
        )
        if len(second_hits) != 1:
            raise ValueError(f"Nonunique mm10-to-mm39 mapping for {name}")
        _, mm39_position, mm39_strand, _ = second_hits[0]
        rows.append({
            "landmark": name,
            "hg38_position_0based": hg38_position,
            "mm39_position_0based": mm39_position,
            "mapped_strand": mm39_strand,
            "distance_to_nearest_mouse_annotated_TSS_bp": min(
                abs(mm39_position - value) for value in mouse_tsses
            ),
            "distance_to_nearest_mouse_CAGE_CTSS_bp": (
                min(abs(mm39_position - value) for value in mouse_ctss)
                if mouse_ctss else ""
            ),
            "mapping_status": "unique_coherent_chain",
        })
    mapped_endpoints = [
        int(row["mm39_position_0based"]) for row in rows
        if row["landmark"] in {
            "promoter_window_start", "promoter_window_end_last_base"
        }
        and row["mapping_status"] == "unique_coherent_chain"
    ]
    mapped_ctss = [
        row for row in rows
        if row["landmark"].startswith("human_CAGE_CTSS_")
        and row["mapping_status"] == "unique_coherent_chain"
    ]
    mapped_annotated = [
        row for row in rows
        if row["landmark"].startswith("human_annotated_TSS:")
        and row["mapping_status"] == "unique_coherent_chain"
    ]
    summary = {
        "human_promoter_window_hg38": {
            "chromosome": "chr2",
            "start_0based": promoter_start,
            "end_0based_exclusive": promoter_end_last_base + 1,
            "strand": "-",
        },
        "mapped_promoter_endpoints": len(mapped_endpoints),
        "mapped_mouse_window_GRCm39": (
            {
                "chromosome": "chr12",
                "start_0based": min(mapped_endpoints),
                "end_0based_exclusive": max(mapped_endpoints) + 1,
                "strand": "+",
                "mouse_curated_full_length_TSSs_in_window": sum(
                    min(mapped_endpoints) <= tss <= max(mapped_endpoints)
                    for tss in mouse_tsses
                ),
            }
            if len(mapped_endpoints) == 2 else None
        ),
        "mapped_human_annotated_TSS_count": len(mapped_annotated),
        "minimum_mapped_human_annotated_TSS_distance_to_mouse_TSS_bp": (
            min(
                int(row["distance_to_nearest_mouse_annotated_TSS_bp"])
                for row in mapped_annotated
            )
            if mapped_annotated else None
        ),
        "mapped_human_CAGE_CTSS_count": len(mapped_ctss),
        "minimum_mapped_human_CAGE_CTSS_distance_to_mouse_TSS_bp": (
            min(
                int(row["distance_to_nearest_mouse_annotated_TSS_bp"])
                for row in mapped_ctss
            )
            if mapped_ctss else None
        ),
        "minimum_mapped_human_CAGE_CTSS_distance_to_mouse_CAGE_CTSS_bp": (
            min(
                int(row["distance_to_nearest_mouse_CAGE_CTSS_bp"])
                for row in mapped_ctss
            )
            if mapped_ctss and mouse_ctss else None
        ),
        "limit": (
            "Syntenic mapping supports locus orthology but does not establish "
            "conserved promoter activity or regulatory effects."
        ),
    }
    return rows, summary


def pair_metric(
    rows: list[dict[str, object]],
    left_species: str,
    right_species: str,
) -> dict[str, object]:
    return next(
        row for row in rows
        if row["region"] == "upstream_2000bp"
        and row["kmer_size"] == 15
        and left_species.replace(" ", "_") in str(row["left_model"])
        and right_species.replace(" ", "_") in str(row["right_model"])
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_targets = read_tsv(TARGETS)
    targets = select_primary_targets(all_targets)
    tss_cluster = select_tss_cluster(all_targets, targets)
    records = promoter_records(targets)
    with (OUT / "canonical_dnmt3a1_promoters.fasta").open("w") as handle:
        for record_id, (_, sequence) in records.items():
            handle.write(f">{record_id}\n{sequence}\n")
    write_tsv(OUT / "canonical_dnmt3a1_targets.tsv", targets)
    write_tsv(OUT / "canonical_dnmt3a1_tss_cluster.tsv", tss_cluster)

    similarity = sequence_similarity(records)
    write_tsv(OUT / "canonical_promoter_sequence_similarity.tsv", similarity)
    nearby, direct = collect_cage(tss_cluster)
    write_tsv(OUT / "canonical_fantom5_nearby_cage_peaks.tsv", nearby)
    write_tsv(OUT / "canonical_direct_tss_evidence.tsv", direct)
    synteny_rows, synteny_summary = human_mouse_synteny(
        tss_cluster, nearby, records
    )
    write_tsv(OUT / "human_mouse_canonical_promoter_synteny.tsv", synteny_rows)

    human_mouse = pair_metric(
        similarity, "Homo sapiens", "Mus musculus"
    )
    tammar_dunnart = pair_metric(
        similarity, "Notamacropus eugenii", "Sminthopsis crassicaudata"
    )
    summary = {
        "representative_promoters": len(targets),
        "canonical_tss_cluster_models": len(tss_cluster),
        "window_definition": (
            "2000 bp upstream plus 500 bp downstream of the selected "
            "full-length-transcript TSS"
        ),
        "representative_rule": (
            "Curated transcript variant 1 for human and mouse; a predicted "
            "full-length outer promoter model for each non-model mammal."
        ),
        "human_mouse_upstream_shared_15mers": human_mouse[
            "shared_unique_kmers"
        ],
        "human_mouse_upstream_longest_exact_block": human_mouse[
            "longest_exact_block_length"
        ],
        "tammar_dunnart_upstream_shared_15mers": tammar_dunnart[
            "shared_unique_kmers"
        ],
        "tammar_dunnart_upstream_longest_exact_block": tammar_dunnart[
            "longest_exact_block_length"
        ],
        "direct_tss_results": direct,
        "human_mouse_synteny": synteny_summary,
        "interpretation_limit": (
            "The analysis tests selected canonical full-length transcript "
            "starts. DNMT3A1 has alternative upstream starts, and exact "
            "sequence similarity does not measure promoter function."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
