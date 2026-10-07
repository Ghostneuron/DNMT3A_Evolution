#!/usr/bin/env python3
"""Test dunnart DNMT3A internal first exons in adult-cerebellum ONT cDNA."""

from __future__ import annotations

import csv
import gzip
import io
import json
import shutil
import subprocess
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path

import ahocorasick
from Bio import Align, SeqIO
from Bio.Seq import Seq

from dunnart_transcriptome_first_exon import junction_seeds, spaced_seeds


ROOT = Path(__file__).resolve().parents[1]
FASTQ = ROOT / "data/raw/marsupial_transcriptomes/ERR15827373.fastq.gz"
GENBANK = (
    ROOT / "data/raw/promoter_tss/"
    "Sminthopsis_crassicaudata_NC_133619.1_DNMT3A.gb"
)
RESULTS = ROOT / "results/06_isoform_evolution"
TARGETS = ("XM_074300067.1", "XM_074300068.1")


@contextmanager
def open_fastq_text(path: Path):
    """Use parallel gzip decompression when pigz is available."""
    pigz = shutil.which("pigz")
    if not pigz:
        with gzip.open(path, "rt") as handle:
            yield handle
        return
    process = subprocess.Popen(
        [pigz, "-dc", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if process.stdout is None:
        raise RuntimeError("pigz did not expose stdout")
    handle = io.TextIOWrapper(process.stdout)
    try:
        yield handle
    finally:
        handle.close()
        stderr = process.stderr.read().decode() if process.stderr else ""
        return_code = process.wait()
        if return_code:
            raise RuntimeError(f"pigz failed ({return_code}): {stderr.strip()}")


def target_exons() -> dict[str, list[str]]:
    record = SeqIO.read(GENBANK, "genbank")
    result = {}
    for target in TARGETS:
        feature = next(
            item for item in record.features
            if item.type == "mRNA"
            and item.qualifiers.get("transcript_id", [""])[0] == target
        )
        result[target] = [
            str(part.extract(record.seq)).upper()
            for part in feature.location.parts
        ]
    return result


def build_seeds(exons: dict[str, list[str]]) -> list[dict[str, str]]:
    rows = []
    for target, parts in exons.items():
        shared_index = 1 if target == TARGETS[0] else 2
        components = [
            ("first_exon", spaced_seeds(parts[0], count=24)),
            ("first_junction", junction_seeds(parts[0], parts[1])),
            (
                "shared_dnmt3a_exon",
                spaced_seeds(parts[shared_index], count=10),
            ),
        ]
        if shared_index == 2:
            components.append(
                ("second_junction", junction_seeds(parts[1], parts[2]))
            )
        for component, component_seeds in components:
            for _, seed in component_seeds:
                for orientation, pattern in (
                    ("forward", seed),
                    ("reverse_complement", str(Seq(seed).reverse_complement())),
                ):
                    rows.append({
                        "target_transcript": target,
                        "component": component,
                        "orientation": orientation,
                        "seed": pattern,
                    })
    return rows


def exact_matches(
    sequence: str,
    automaton: ahocorasick.Automaton,
    seed_lookup: dict[str, list[dict[str, str]]],
) -> dict[tuple[str, str, str], int]:
    matched_seed_strings = {seed for _, seed in automaton.iter(sequence.upper())}
    matches: dict[tuple[str, str, str], int] = defaultdict(int)
    for seed in matched_seed_strings:
        for row in seed_lookup[seed]:
            matches[(
                row["target_transcript"],
                row["component"],
                row["orientation"],
            )] += 1
    return dict(matches)


def main() -> None:
    if not FASTQ.exists():
        raise SystemExit(f"Missing FASTQ: {FASTQ}")
    exons = target_exons()
    seeds = build_seeds(exons)
    seed_lookup: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in seeds:
        seed_lookup[row["seed"]].append(row)
    automaton = ahocorasick.Automaton()
    for seed in seed_lookup:
        automaton.add_word(seed, seed)
    automaton.make_automaton()

    evidence_rows = []
    matched_reads: list[tuple[str, str, str]] = []
    total_reads = 0
    orphan_header_ids = []
    with open_fastq_text(FASTQ) as handle:
        pending_header = ""
        while True:
            header = pending_header or handle.readline()
            pending_header = ""
            if not header:
                break
            if not header.startswith("@"):
                raise ValueError(
                    f"Expected FASTQ header after {total_reads} reads: "
                    f"{header[:80]!r}"
                )
            sequence_line = handle.readline()
            if sequence_line.startswith("@"):
                # The submitted source contains at least one header-only
                # record. Preserve an audit and resume at the following header.
                orphan_header_ids.append(header[1:].split()[0])
                pending_header = sequence_line
                continue
            sequence = sequence_line.strip().upper()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not plus.startswith("+") or len(quality) != len(sequence):
                raise ValueError(
                    f"Malformed FASTQ record {header[1:].split()[0]} "
                    f"after {total_reads} complete reads"
                )
            title = header[1:].rstrip()
            sequence = sequence.upper()
            total_reads += 1
            matches = exact_matches(sequence, automaton, seed_lookup)
            if not matches:
                if total_reads % 1_000_000 == 0:
                    print(f"scanned_reads={total_reads}", flush=True)
                continue
            read_id = title.split()[0]
            matched_reads.append((f"@{title}", sequence, quality))
            for (target, component, orientation), count in matches.items():
                evidence_rows.append({
                    "run_accession": "ERR15827373",
                    "target_transcript": target,
                    "read_id": read_id,
                    "read_length": len(sequence),
                    "component": component,
                    "orientation": orientation,
                    "distinct_exact_35mer_hits": count,
                })
            if total_reads % 1_000_000 == 0:
                print(f"scanned_reads={total_reads}", flush=True)

    evidence_rows.sort(key=lambda row: (
        row["target_transcript"], row["read_id"], row["component"],
        row["orientation"],
    ))
    evidence_output = RESULTS / "dunnart_cerebellum_ont_first_exon_evidence.tsv"
    with evidence_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(evidence_rows[0])
        )
        writer.writeheader()
        writer.writerows(evidence_rows)
    with (RESULTS / "dunnart_cerebellum_ont_dnmt3a_seed_reads.fastq").open(
        "w"
    ) as handle:
        for header, sequence, quality in matched_reads:
            handle.write(f"{header}\n{sequence}\n+\n{quality}\n")

    read_sequences = {
        header[1:].split()[0]: sequence
        for header, sequence, _ in matched_reads
    }
    components: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in evidence_rows:
        components[(row["target_transcript"], row["read_id"])].add(
            row["component"]
        )

    aligner = Align.PairwiseAligner(
        mode="local",
        match_score=2,
        mismatch_score=-1,
        open_gap_score=-5,
        extend_gap_score=-1,
    )
    alignment_rows = []
    target_summary = {}
    for target in TARGETS:
        first_junction_reads = sorted(
            read_id for (row_target, read_id), values in components.items()
            if row_target == target and "first_junction" in values
        )
        first_and_shared_reads = sorted(
            read_id for (row_target, read_id), values in components.items()
            if row_target == target
            and {"first_exon", "shared_dnmt3a_exon"} <= values
        )
        model_sequence = "".join(exons[target])
        for read_id in first_junction_reads:
            alternatives = []
            for orientation, oriented_read in (
                ("forward", read_sequences[read_id]),
                (
                    "reverse_complement",
                    str(Seq(read_sequences[read_id]).reverse_complement()),
                ),
            ):
                alignment = aligner.align(model_sequence, oriented_read)[0]
                alternatives.append((alignment.score, orientation, alignment))
            score, orientation, alignment = max(
                alternatives, key=lambda item: item[0]
            )
            counts = alignment.counts()
            coordinates = alignment.coordinates
            alignment_rows.append({
                "target_transcript": target,
                "read_id": read_id,
                "read_length": len(read_sequences[read_id]),
                "orientation_to_model": orientation,
                "alignment_score": f"{score:.1f}",
                "model_start_0based": int(coordinates[0, 0]),
                "model_end_0based_exclusive": int(coordinates[0, -1]),
                "read_start_0based": int(coordinates[1, 0]),
                "read_end_0based_exclusive": int(coordinates[1, -1]),
                "aligned_identical_bases": counts.identities,
                "aligned_nongap_bases": counts.aligned,
                "aligned_identity": f"{counts.identities / counts.aligned:.4f}",
                "reaches_annotated_transcript_base_1": str(
                    int(coordinates[0, 0]) == 0
                ).lower(),
            })
        target_summary[target] = {
            "reads_with_exact_first_junction_seed": first_junction_reads,
            "reads_with_first_and_shared_exon_seeds": first_and_shared_reads,
        }

    alignment_output = (
        RESULTS / "dunnart_cerebellum_ont_junction_read_alignments.tsv"
    )
    with alignment_output.open("w", newline="") as handle:
        fieldnames = [
            "target_transcript", "read_id", "read_length",
            "orientation_to_model", "alignment_score", "model_start_0based",
            "model_end_0based_exclusive", "read_start_0based",
            "read_end_0based_exclusive", "aligned_identical_bases",
            "aligned_nongap_bases", "aligned_identity",
            "reaches_annotated_transcript_base_1",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(alignment_rows)

    for target in TARGETS:
        target_summary[target][
            "junction_reads_reaching_annotated_transcript_base_1"
        ] = [
            row["read_id"] for row in alignment_rows
            if row["target_transcript"] == target
            and row["reaches_annotated_transcript_base_1"] == "true"
        ]
    summary = {
        "run_accession": "ERR15827373",
        "sample": "adult female cerebellum",
        "platform": "Oxford Nanopore PromethION cDNA",
        "fastq": str(FASTQ.relative_to(ROOT)),
        "fastq_expected_md5": "8866eda56108ec2620b7b2c341cc3b76",
        "reads_scanned": total_reads,
        "header_only_records_skipped": len(orphan_header_ids),
        "header_only_record_ids": orphan_header_ids,
        "reads_with_any_exact_35mer_seed": len(matched_reads),
        "targets": target_summary,
        "interpretation_limit": (
            "Junction-spanning cerebellum cDNA establishes brain-context "
            "expression of an internal transcript, but non-cap-enriched read "
            "starts do not establish the exact TSS."
        ),
    }
    (RESULTS / "dunnart_cerebellum_ont_first_exon_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
