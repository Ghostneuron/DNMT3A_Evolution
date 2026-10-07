#!/usr/bin/env python3
"""Build a competitive DNMT3A target/background reference for Salmon."""

from __future__ import annotations

import csv
from pathlib import Path

import ahocorasick
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parents[1]
GENBANK = (
    ROOT / "data/raw/promoter_tss/"
    "Sminthopsis_crassicaudata_NC_133619.1_DNMT3A.gb"
)
TRINITY = (
    ROOT / "data/raw/marsupial_transcriptomes/"
    "Dunnart_trinity_k29_output.Trinity.fasta"
)
REFERENCE_DIR = ROOT / "data/processed/dunnart_salmon_reference"
TARGETS = REFERENCE_DIR / "dnmt3a_refseq_targets.fasta"
BACKGROUND = REFERENCE_DIR / "trinity_non_dnmt3a_decoys.fasta"
DECOYS = REFERENCE_DIR / "decoys.txt"
TARGET_AUDIT = REFERENCE_DIR / "dnmt3a_target_audit.tsv"
EXCLUDED_AUDIT = REFERENCE_DIR / "excluded_dnmt3a_like_trinity.tsv"
KMER_LENGTH = 31
EXCLUSION_THRESHOLD = 100


def kmers(sequence: str, length: int = KMER_LENGTH) -> set[str]:
    sequence = sequence.upper()
    return {
        sequence[index:index + length]
        for index in range(len(sequence) - length + 1)
        if set(sequence[index:index + length]) <= {"A", "C", "G", "T"}
    }


def shared_kmer_count(
    sequence: str,
    target_kmers: set[str],
    length: int = KMER_LENGTH,
) -> int:
    return len(kmers(sequence, length=length) & target_kmers)


def kmer_automaton(target_kmers: set[str]) -> ahocorasick.Automaton:
    automaton = ahocorasick.Automaton()
    for kmer in target_kmers:
        automaton.add_word(kmer, kmer)
    automaton.make_automaton()
    return automaton


def automaton_shared_kmer_count(
    sequence: str,
    automaton: ahocorasick.Automaton,
) -> int:
    return len({kmer for _, kmer in automaton.iter(sequence.upper())})


def target_records() -> list[SeqRecord]:
    record = SeqIO.read(GENBANK, "genbank")
    targets = []
    for feature in record.features:
        if feature.type != "mRNA":
            continue
        transcript_id = feature.qualifiers.get("transcript_id", [""])[0]
        targets.append(SeqRecord(
            feature.extract(record.seq),
            id=transcript_id,
            description="DNMT3A RefSeq mRNA model",
        ))
    return sorted(targets, key=lambda item: item.id)


def main() -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    targets = target_records()
    SeqIO.write(targets, TARGETS, "fasta")
    target_kmers = set().union(*(kmers(str(record.seq)) for record in targets))
    automaton = kmer_automaton(target_kmers)

    with TARGET_AUDIT.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["transcript_id", "length", "unique_31mers"])
        for record in targets:
            writer.writerow([record.id, len(record.seq), len(kmers(str(record.seq)))])

    excluded = []
    decoy_count = 0
    with BACKGROUND.open("w") as fasta_handle, DECOYS.open("w") as decoy_handle:
        for record in SeqIO.parse(TRINITY, "fasta"):
            count = automaton_shared_kmer_count(str(record.seq), automaton)
            if count >= EXCLUSION_THRESHOLD:
                excluded.append((record.id, len(record.seq), count))
                continue
            decoy_id = f"DECOY|{record.id}"
            record.id = decoy_id
            record.name = decoy_id
            record.description = ""
            SeqIO.write(record, fasta_handle, "fasta")
            decoy_handle.write(decoy_id + "\n")
            decoy_count += 1

    with EXCLUDED_AUDIT.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "trinity_transcript_id", "length", "distinct_target_31mer_hits",
        ])
        writer.writerows(sorted(excluded, key=lambda row: (-row[2], row[0])))

    print(f"DNMT3A targets: {len(targets)}")
    print(f"Background decoys: {decoy_count}")
    print(f"Excluded DNMT3A-like Trinity transcripts: {len(excluded)}")


if __name__ == "__main__":
    main()
