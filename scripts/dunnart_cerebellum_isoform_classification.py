#!/usr/bin/env python3
"""Classify retained dunnart cerebellum reads by paralog and DNMT3A isoform."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Seq import Seq

from dunnart_transcriptome_first_exon import junction_seeds


ROOT = Path(__file__).resolve().parents[1]
READS = (
    ROOT / "results/06_isoform_evolution/"
    "dunnart_cerebellum_ont_dnmt3a_seed_reads.fastq"
)
GENBANK = (
    ROOT / "data/raw/promoter_tss/"
    "Sminthopsis_crassicaudata_NC_133619.1_DNMT3A.gb"
)
DNMT3A_PROTEINS = ROOT / "data/raw/ncbi_dataset/data/protein.faa"
DNMT3B_PROTEINS = (
    ROOT / "data/raw/marsupial_transcriptomes/dunnart_DNMT3B_proteins.fasta"
)
RESULTS = ROOT / "results/06_isoform_evolution"
FULL_LENGTH = "XM_074300059.1"
INTERNAL_MODELS = ("XM_074300067.1", "XM_074300068.1")


def transcript_exons(record, transcript_id: str) -> list[str]:
    feature = next(
        item for item in record.features
        if item.type == "mRNA"
        and item.qualifiers.get("transcript_id", [""])[0] == transcript_id
    )
    return [
        str(part.extract(record.seq)).upper()
        for part in feature.location.parts
    ]


def both_orientation_junction_seeds(parts: list[str]) -> set[str]:
    forward = {seed for _, seed in junction_seeds(parts[0], parts[1])}
    return forward | {str(Seq(seed).reverse_complement()) for seed in forward}


def exact_seed_count(sequence: str, seeds: set[str]) -> int:
    sequence = sequence.upper()
    return sum(seed in sequence for seed in seeds)


def classify_isoform(
    full_length_hits: int,
    internal_067_hits: int,
    internal_068_hits: int,
) -> str:
    if internal_067_hits:
        return "internal_XM_074300067"
    if internal_068_hits:
        return "internal_XM_074300068"
    if full_length_hits:
        return "full_length_XM_074300059"
    return "downstream_DNMT3A_only"


def translated_orfs(sequence: str, minimum: int = 40):
    for orientation, nucleotide in (
        ("forward", sequence),
        ("reverse_complement", str(Seq(sequence).reverse_complement())),
    ):
        for frame in range(3):
            usable = nucleotide[frame:]
            usable = usable[:len(usable) - len(usable) % 3]
            for fragment in str(Seq(usable).translate()).split("*"):
                if len(fragment) >= minimum:
                    yield orientation, frame, fragment


def best_protein_score(aligner, orfs, reference: str):
    values = []
    for orientation, frame, amino_acids in orfs:
        alignment = aligner.align(amino_acids, reference)[0]
        counts = alignment.counts()
        values.append((
            alignment.score,
            counts.identities,
            counts.aligned,
            orientation,
            frame,
            len(amino_acids),
        ))
    return max(values) if values else (0, 0, 0, "", 0, 0)


def main() -> None:
    record = SeqIO.read(GENBANK, "genbank")
    junction_sets = {
        FULL_LENGTH: both_orientation_junction_seeds(
            transcript_exons(record, FULL_LENGTH)
        ),
        **{
            target: both_orientation_junction_seeds(
                transcript_exons(record, target)
            )
            for target in INTERNAL_MODELS
        },
    }
    dnmt3a = next(
        str(item.seq) for item in SeqIO.parse(DNMT3A_PROTEINS, "fasta")
        if item.id == "XP_074156160.1"
    )
    dnmt3b = next(
        str(item.seq) for item in SeqIO.parse(DNMT3B_PROTEINS, "fasta")
        if item.id == "XP_074148839.1"
    )
    aligner = Align.PairwiseAligner(
        mode="local",
        match_score=2,
        mismatch_score=-1,
        open_gap_score=-5,
        extend_gap_score=-1,
    )
    rows = []
    for read in SeqIO.parse(READS, "fastq"):
        sequence = str(read.seq)
        orfs = list(translated_orfs(sequence))
        score_a = best_protein_score(aligner, orfs, dnmt3a)
        score_b = best_protein_score(aligner, orfs, dnmt3b)
        if score_a[0] - score_b[0] >= 20:
            paralog = "DNMT3A"
        elif score_b[0] - score_a[0] >= 20:
            paralog = "DNMT3B"
        else:
            paralog = "ambiguous"
        full_hits = exact_seed_count(sequence, junction_sets[FULL_LENGTH])
        internal_067_hits = exact_seed_count(
            sequence, junction_sets[INTERNAL_MODELS[0]]
        )
        internal_068_hits = exact_seed_count(
            sequence, junction_sets[INTERNAL_MODELS[1]]
        )
        rows.append({
            "read_id": read.id,
            "read_length": len(read.seq),
            "paralog_classification": paralog,
            "best_DNMT3A_protein_score": f"{score_a[0]:.1f}",
            "best_DNMT3B_protein_score": f"{score_b[0]:.1f}",
            "protein_score_margin_A_minus_B": f"{score_a[0] - score_b[0]:.1f}",
            "full_length_first_junction_exact_35mer_hits": full_hits,
            "internal_067_first_junction_exact_35mer_hits": internal_067_hits,
            "internal_068_first_junction_exact_35mer_hits": internal_068_hits,
            "isoform_evidence_classification": classify_isoform(
                full_hits, internal_067_hits, internal_068_hits
            ),
        })
    output = RESULTS / "dunnart_cerebellum_ont_isoform_classification.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(rows[0])
        )
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "candidate_reads": len(rows),
        "paralog_classification_counts": dict(
            Counter(row["paralog_classification"] for row in rows)
        ),
        "isoform_evidence_counts": dict(
            Counter(row["isoform_evidence_classification"] for row in rows)
        ),
        "full_length_first_junction_reads": sum(
            bool(row["full_length_first_junction_exact_35mer_hits"])
            for row in rows
        ),
        "internal_first_junction_reads": sum(
            bool(row["internal_067_first_junction_exact_35mer_hits"])
            or bool(row["internal_068_first_junction_exact_35mer_hits"])
            for row in rows
        ),
        "inference": (
            "DNMT3A is expressed in adult dunnart cerebellum, with direct "
            "support for a full-length-model first junction but no exact "
            "support for either predicted internal first junction in this run."
        ),
        "limit": (
            "Failure to observe an internal junction in one adult female "
            "cerebellum replicate is not evidence of organism-wide absence."
        ),
    }
    (RESULTS / "dunnart_cerebellum_ont_isoform_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
