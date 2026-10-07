#!/usr/bin/env python3
"""Descriptive exact-sequence similarity among locus-anchored promoter windows."""

from __future__ import annotations

import csv
import json
from difflib import SequenceMatcher
from itertools import combinations
from pathlib import Path

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "results/06_isoform_evolution/primary_internal_start_promoters.fasta"
RESULTS = ROOT / "results/06_isoform_evolution"
TSS_OFFSET = 2000
REGIONS = {
    "upstream_2000bp": (0, 2000),
    "downstream_500bp": (2000, 2500),
}
KMER_SIZES = (11, 15, 19)


def kmers(sequence: str, size: int) -> set[str]:
    return {
        sequence[index:index + size]
        for index in range(len(sequence) - size + 1)
        if "N" not in sequence[index:index + size]
    }


def main() -> None:
    records = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(FASTA, "fasta")
    }
    rows = []
    for (left_id, left), (right_id, right) in combinations(records.items(), 2):
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
                        "Descriptive locus-constrained similarity; repeats and "
                        "short motifs can match without regulatory conservation"
                    ),
                })

    with (RESULTS / "promoter_exact_sequence_similarity.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    upstream_15 = [
        row for row in rows
        if row["region"] == "upstream_2000bp" and row["kmer_size"] == 15
    ]
    def pair(left_prefix: str, right_prefix: str):
        return next(
            row for row in upstream_15
            if (
                row["left_model"].startswith(left_prefix)
                and row["right_model"].startswith(right_prefix)
            ) or (
                row["left_model"].startswith(right_prefix)
                and row["right_model"].startswith(left_prefix)
            )
        )

    human_mouse = pair("Homo_sapiens", "Mus_musculus")
    tammar_dunnart = [
        row for row in upstream_15
        if (
            "Notamacropus_eugenii" in row["left_model"]
            and "Sminthopsis_crassicaudata" in row["right_model"]
        ) or (
            "Notamacropus_eugenii" in row["right_model"]
            and "Sminthopsis_crassicaudata" in row["left_model"]
        )
    ]
    human_marsupial = [
        row for row in upstream_15
        if (
            row["left_model"].startswith("Homo_sapiens")
            and (
                row["right_model"].startswith("Notamacropus_eugenii")
                or row["right_model"].startswith("Sminthopsis_crassicaudata")
            )
        )
    ]
    summary = {
        "window_definition": "2000 bp upstream plus 500 bp downstream of annotated TSS",
        "human_mouse_upstream_shared_15mers": human_mouse["shared_unique_kmers"],
        "human_mouse_upstream_longest_exact_block": human_mouse[
            "longest_exact_block_length"
        ],
        "tammar_dunnart_upstream_shared_15mer_range": [
            min(row["shared_unique_kmers"] for row in tammar_dunnart),
            max(row["shared_unique_kmers"] for row in tammar_dunnart),
        ],
        "tammar_dunnart_upstream_longest_exact_block_range": [
            min(row["longest_exact_block_length"] for row in tammar_dunnart),
            max(row["longest_exact_block_length"] for row in tammar_dunnart),
        ],
        "human_marsupial_upstream_shared_15mer_range": [
            min(row["shared_unique_kmers"] for row in human_marsupial),
            max(row["shared_unique_kmers"] for row in human_marsupial),
        ],
        "inference": (
            "Upstream sequence conservation is clear within human-mouse and "
            "especially between the two marsupial loci, but weak across the "
            "human-marsupial comparison at exact-kmer resolution."
        ),
        "limit": (
            "This is not a genome-wide orthology map or a promoter-function "
            "test; downstream similarity is partly driven by transcribed exons."
        ),
    }
    (RESULTS / "promoter_exact_sequence_similarity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
