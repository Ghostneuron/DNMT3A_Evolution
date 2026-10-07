#!/usr/bin/env python3
"""Search the public dunnart Trinity assembly for DNMT3A internal first exons."""

from __future__ import annotations

import csv
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[1]
TRANSCRIPTOME = (
    ROOT / "data/raw/marsupial_transcriptomes/"
    "Dunnart_trinity_k29_output.Trinity.fasta"
)
GENBANK = (
    ROOT / "data/raw/promoter_tss/"
    "Sminthopsis_crassicaudata_NC_133619.1_DNMT3A.gb"
)
RESULTS = ROOT / "results/06_isoform_evolution"
TARGETS = ("XM_074300067.1", "XM_074300068.1")
SEED_LENGTH = 35


def sequence_entropy(sequence: str) -> float:
    length = len(sequence)
    return -sum(
        (sequence.count(base) / length) * math.log2(sequence.count(base) / length)
        for base in "ACGT"
        if sequence.count(base)
    )


def spaced_seeds(sequence: str, count: int = 16) -> list[tuple[int, str]]:
    sequence = sequence.upper()
    maximum = len(sequence) - SEED_LENGTH
    if maximum < 0:
        return []
    starts = sorted({
        round(index * maximum / max(count - 1, 1))
        for index in range(count)
    })
    return [
        (start, sequence[start:start + SEED_LENGTH])
        for start in starts
        if set(sequence[start:start + SEED_LENGTH]) <= set("ACGT")
        and sequence_entropy(sequence[start:start + SEED_LENGTH]) >= 1.55
        and max(sequence[start:start + SEED_LENGTH].count(base) for base in "ACGT")
        / SEED_LENGTH < 0.60
    ]


def junction_seeds(left: str, right: str) -> list[tuple[int, str]]:
    combined = left + right
    boundary = len(left)
    rows = []
    for left_bases in (10, 14, 17, 21, 25):
        start = boundary - left_bases
        seed = combined[start:start + SEED_LENGTH].upper()
        if len(seed) == SEED_LENGTH and sequence_entropy(seed) >= 1.55:
            rows.append((start, seed))
    return rows


def target_exons() -> dict[str, list[str]]:
    record = SeqIO.read(GENBANK, "genbank")
    exons: dict[str, list[str]] = {}
    for target in TARGETS:
        feature = next(
            item for item in record.features
            if item.type == "mRNA"
            and item.qualifiers.get("transcript_id", [""])[0] == target
        )
        exons[target] = [
            str(part.extract(record.seq)).upper()
            for part in feature.location.parts
        ]
    return exons


def build_seed_manifest() -> list[dict[str, str | int]]:
    exons = target_exons()
    rows: list[dict[str, str | int]] = []
    for target, parts in exons.items():
        shared_index = 1 if target == TARGETS[0] else 2
        components = [
            ("first_exon", spaced_seeds(parts[0])),
            ("first_junction", junction_seeds(parts[0], parts[1])),
            ("shared_dnmt3a_exon", spaced_seeds(parts[shared_index], count=8)),
        ]
        if shared_index == 2:
            components.append(
                ("second_junction", junction_seeds(parts[1], parts[2]))
            )
        for component, seeds in components:
            for query_offset, seed in seeds:
                for orientation, pattern in (
                    ("forward", seed),
                    ("reverse_complement", str(Seq(seed).reverse_complement())),
                ):
                    rows.append({
                        "target_transcript": target,
                        "component": component,
                        "query_offset": query_offset,
                        "orientation": orientation,
                        "seed": pattern,
                    })
    return rows


def run_ripgrep(seed_file: Path) -> list[tuple[int, str]]:
    process = subprocess.run(
        [
            "rg", "--line-number", "--only-matching", "--fixed-strings",
            "--file", str(seed_file), str(TRANSCRIPTOME),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if process.returncode not in (0, 1):
        raise RuntimeError(process.stderr.strip())
    hits = []
    for line in process.stdout.splitlines():
        line_number, matched = line.split(":", 1)
        hits.append((int(line_number), matched.upper()))
    return hits


def fasta_headers_for_lines(line_numbers: set[int]) -> dict[int, str]:
    matches: dict[int, str] = {}
    header = ""
    with TRANSCRIPTOME.open() as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                header = line[1:].strip()
            if number in line_numbers:
                matches[number] = header
    return matches


def fasta_records_by_id(record_ids: set[str]) -> dict[str, Seq]:
    records: dict[str, Seq] = {}
    for record in SeqIO.parse(TRANSCRIPTOME, "fasta"):
        if record.id in record_ids:
            records[record.id] = record.seq
            if len(records) == len(record_ids):
                break
    return records


def main() -> None:
    if not TRANSCRIPTOME.exists():
        raise SystemExit(f"Missing transcriptome: {TRANSCRIPTOME}")
    RESULTS.mkdir(parents=True, exist_ok=True)
    seed_rows = build_seed_manifest()
    seed_manifest = RESULTS / "dunnart_first_exon_seed_manifest.tsv"
    with seed_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(seed_rows[0])
        )
        writer.writeheader()
        writer.writerows(seed_rows)

    seed_file = RESULTS / "dunnart_first_exon_search_seeds.txt"
    seed_file.write_text(
        "\n".join(sorted({str(row["seed"]) for row in seed_rows})) + "\n"
    )
    hits = run_ripgrep(seed_file)
    headers = fasta_headers_for_lines({line for line, _ in hits})
    seed_lookup: dict[str, list[dict[str, str | int]]] = defaultdict(list)
    for row in seed_rows:
        seed_lookup[str(row["seed"])].append(row)

    evidence: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    for line_number, matched in hits:
        trinity_id = headers[line_number].split()[0]
        for seed_row in seed_lookup[matched]:
            key = (
                str(seed_row["target_transcript"]),
                trinity_id,
                str(seed_row["component"]),
                str(seed_row["orientation"]),
            )
            evidence[key].add(matched)

    evidence_rows = [
        {
            "target_transcript": target,
            "trinity_transcript": trinity,
            "component": component,
            "orientation": orientation,
            "distinct_exact_35mer_hits": len(seeds),
            "evidence_interpretation": (
                "splice-junction support"
                if "junction" in component
                else "exonic sequence support"
            ),
        }
        for (target, trinity, component, orientation), seeds in evidence.items()
    ]
    evidence_rows.sort(key=lambda row: (
        row["target_transcript"],
        row["trinity_transcript"],
        row["component"],
        row["orientation"],
    ))

    # ripgrep identifies candidate records quickly, but a biological 35-mer can
    # cross a FASTA display-line boundary. Re-evaluate every seed against each
    # complete candidate sequence before classifying exon-junction evidence.
    candidate_records = fasta_records_by_id({
        str(row["trinity_transcript"]) for row in evidence_rows
    })
    evidence = defaultdict(set)
    for trinity_id, sequence in candidate_records.items():
        complete_sequence = str(sequence).upper()
        for seed_row in seed_rows:
            seed = str(seed_row["seed"])
            if seed not in complete_sequence:
                continue
            key = (
                str(seed_row["target_transcript"]),
                trinity_id,
                str(seed_row["component"]),
                str(seed_row["orientation"]),
            )
            evidence[key].add(seed)
    evidence_rows = [
        {
            "target_transcript": target,
            "trinity_transcript": trinity,
            "component": component,
            "orientation": orientation,
            "distinct_exact_35mer_hits": len(seeds),
            "evidence_interpretation": (
                "splice-junction support"
                if "junction" in component
                else "exonic sequence support"
            ),
        }
        for (target, trinity, component, orientation), seeds in evidence.items()
    ]
    evidence_rows.sort(key=lambda row: (
        row["target_transcript"],
        row["trinity_transcript"],
        row["component"],
        row["orientation"],
    ))
    output = RESULTS / "dunnart_transcriptome_first_exon_evidence.tsv"
    with output.open("w", newline="") as handle:
        fieldnames = [
            "target_transcript", "trinity_transcript", "component",
            "orientation", "distinct_exact_35mer_hits",
            "evidence_interpretation",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(evidence_rows)

    by_target: dict[str, dict[str, set[str]]] = {
        target: defaultdict(set) for target in TARGETS
    }
    for row in evidence_rows:
        by_target[str(row["target_transcript"])][
            str(row["trinity_transcript"])
        ].add(str(row["component"]))
    summary_targets = {}
    for target, transcripts in by_target.items():
        junction_supported = sorted(
            transcript for transcript, components in transcripts.items()
            if any("junction" in component for component in components)
        )
        first_and_shared = sorted(
            transcript for transcript, components in transcripts.items()
            if {"first_exon", "shared_dnmt3a_exon"} <= components
        )
        summary_targets[target] = {
            "trinity_transcripts_with_any_seed": len(transcripts),
            "trinity_transcripts_with_exact_junction_seed": junction_supported,
            "trinity_transcripts_with_first_and_shared_exon_seeds": first_and_shared,
        }

    selected_ids = {
        transcript
        for target_summary in summary_targets.values()
        for key in (
            "trinity_transcripts_with_exact_junction_seed",
            "trinity_transcripts_with_first_and_shared_exon_seeds",
        )
        for transcript in target_summary[key]
    }
    selected_records = fasta_records_by_id(selected_ids)
    candidate_fasta = RESULTS / "dunnart_dnmt3a_candidate_transcripts.fasta"
    with candidate_fasta.open("w") as handle:
        for record_id in sorted(selected_records):
            handle.write(f">{record_id}\n")
            sequence = str(selected_records[record_id])
            handle.write(
                "\n".join(
                    sequence[start:start + 80]
                    for start in range(0, len(sequence), 80)
                ) + "\n"
            )
    summary = {
        "transcriptome": str(TRANSCRIPTOME.relative_to(ROOT)),
        "transcriptome_figshare_file_id": 42751783,
        "transcriptome_expected_md5": "5defd5fb1a32057b42df15526e31c847",
        "seed_length": SEED_LENGTH,
        "candidate_transcript_fasta": str(candidate_fasta.relative_to(ROOT)),
        "candidate_transcripts_written": len(selected_records),
        "targets": summary_targets,
        "interpretation_limit": (
            "An assembled cDNA match supports transcript structure but cannot "
            "identify the exact capped transcription-start nucleotide."
        ),
    }
    (RESULTS / "dunnart_transcriptome_first_exon_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
