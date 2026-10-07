#!/usr/bin/env python3
"""Scan paired dunnart neocortex RNA-seq for DNMT3A isoform junctions."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import ahocorasick
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from dunnart_cerebellum_ont_first_exon import open_fastq_text
from dunnart_transcriptome_first_exon import junction_seeds


ROOT = Path(__file__).resolve().parents[1]
GENBANK = (
    ROOT / "data/raw/promoter_tss/"
    "Sminthopsis_crassicaudata_NC_133619.1_DNMT3A.gb"
)
RESULTS = ROOT / "results/06_isoform_evolution"
FULL_LENGTH = "XM_074300059.1"
INTERNAL_067 = "XM_074300067.1"
INTERNAL_068 = "XM_074300068.1"


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


def build_seed_rows() -> list[dict[str, str]]:
    record = SeqIO.read(GENBANK, "genbank")
    exons = {
        target: transcript_exons(record, target)
        for target in (FULL_LENGTH, INTERNAL_067, INTERNAL_068)
    }
    components = [
        (
            "full_length_first_junction",
            FULL_LENGTH,
            exons[FULL_LENGTH][0],
            exons[FULL_LENGTH][1],
        ),
        (
            "internal_067_first_junction",
            INTERNAL_067,
            exons[INTERNAL_067][0],
            exons[INTERNAL_067][1],
        ),
        (
            "internal_068_first_junction",
            INTERNAL_068,
            exons[INTERNAL_068][0],
            exons[INTERNAL_068][1],
        ),
        (
            "common_downstream_junction",
            INTERNAL_067,
            exons[INTERNAL_067][1],
            exons[INTERNAL_067][2],
        ),
    ]
    rows = []
    for component, model, left, right in components:
        for _, seed in junction_seeds(left, right):
            for orientation, pattern in (
                ("forward", seed),
                ("reverse_complement", str(Seq(seed).reverse_complement())),
            ):
                rows.append({
                    "component": component,
                    "model_transcript": model,
                    "orientation": orientation,
                    "seed": pattern,
                })
    return rows


def exact_matches(
    sequence: str,
    automaton: ahocorasick.Automaton,
    lookup: dict[str, list[dict[str, str]]],
) -> dict[tuple[str, str, str], int]:
    matched = {seed for _, seed in automaton.iter(sequence.upper())}
    counts: dict[tuple[str, str, str], int] = defaultdict(int)
    for seed in matched:
        for row in lookup[seed]:
            counts[(
                row["component"],
                row["model_transcript"],
                row["orientation"],
            )] += 1
    return dict(counts)


def normalized_pair_id(title: str) -> str:
    read_id = title.split()[0]
    if read_id.endswith("/1") or read_id.endswith("/2"):
        return read_id[:-2]
    return read_id


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="SRR13036046")
    parser.add_argument("--sample", default="P20_neocortex_rep1")
    parser.add_argument(
        "--sample-description",
        default="P20 neocortex replicate 1, specimen 666",
    )
    parser.add_argument("--mate1", type=Path)
    parser.add_argument("--mate2", type=Path)
    parser.add_argument(
        "--output-prefix",
        default="dunnart_P20_neocortex_rep1_junction",
    )
    return parser.parse_args()


def transcript_seed_specificity(
    record, seed_rows: list[dict[str, str]]
) -> list[dict[str, object]]:
    """Report which annotated DNMT3A transcript models contain each seed."""
    transcripts = {}
    for feature in record.features:
        if feature.type != "mRNA":
            continue
        transcript_id = feature.qualifiers.get("transcript_id", [""])[0]
        transcripts[transcript_id] = "".join(
            str(part.extract(record.seq)).upper()
            for part in feature.location.parts
        )
    rows = []
    for row in seed_rows:
        if row["orientation"] != "forward":
            continue
        seed = row["seed"]
        hits = sorted(
            transcript_id
            for transcript_id, sequence in transcripts.items()
            if seed in sequence
        )
        rows.append({
            "component": row["component"],
            "model_transcript": row["model_transcript"],
            "seed": seed,
            "matching_annotated_transcripts": ",".join(hits),
            "matching_transcript_count": len(hits),
        })
    return rows


def main() -> None:
    args = parse_args()
    fastqs = (
        args.mate1 or ROOT / f"data/raw/marsupial_transcriptomes/{args.run}_1.fastq.gz",
        args.mate2 or ROOT / f"data/raw/marsupial_transcriptomes/{args.run}_2.fastq.gz",
    )
    missing = [str(path) for path in fastqs if not path.exists()]
    if missing:
        raise SystemExit(f"Missing FASTQ files: {missing}")
    seed_rows = build_seed_rows()
    lookup: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in seed_rows:
        lookup[row["seed"]].append(row)
    automaton = ahocorasick.Automaton()
    for seed in lookup:
        automaton.add_word(seed, seed)
    automaton.make_automaton()

    evidence_rows = []
    read_counts = {}
    for mate, path in enumerate(fastqs, start=1):
        count = 0
        with open_fastq_text(path) as handle:
            for title, sequence, _ in FastqGeneralIterator(handle):
                count += 1
                matches = exact_matches(sequence, automaton, lookup)
                if not matches:
                    continue
                for (component, model, orientation), seed_count in matches.items():
                    evidence_rows.append({
                        "run_accession": args.run,
                        "sample": args.sample,
                        "mate": mate,
                        "read_id": title.split()[0],
                        "pair_id": normalized_pair_id(title),
                        "read_length": len(sequence),
                        "component": component,
                        "model_transcript": model,
                        "orientation": orientation,
                        "distinct_exact_35mer_hits": seed_count,
                    })
        read_counts[f"mate_{mate}"] = count

    evidence_rows.sort(key=lambda row: (
        row["component"], row["pair_id"], int(row["mate"])
    ))
    output = RESULTS / f"{args.output_prefix}_evidence.tsv"
    with output.open("w", newline="") as handle:
        fieldnames = [
            "run_accession", "sample", "mate", "read_id", "pair_id",
            "read_length", "component", "model_transcript", "orientation",
            "distinct_exact_35mer_hits",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(evidence_rows)

    pairs_by_component: dict[str, set[str]] = defaultdict(set)
    for row in evidence_rows:
        pairs_by_component[row["component"]].add(row["pair_id"])
    components = (
        "full_length_first_junction",
        "internal_067_first_junction",
        "internal_068_first_junction",
        "common_downstream_junction",
    )
    record = SeqIO.read(GENBANK, "genbank")
    specificity_rows = transcript_seed_specificity(record, seed_rows)
    specificity_output = RESULTS / f"{args.output_prefix}_seed_specificity.tsv"
    with specificity_output.open("w", newline="") as handle:
        fieldnames = [
            "component", "model_transcript", "seed",
            "matching_annotated_transcripts", "matching_transcript_count",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(specificity_rows)

    matching_models = defaultdict(set)
    for row in specificity_rows:
        matching_models[row["component"]].update(
            str(row["matching_annotated_transcripts"]).split(",")
        )
    summary = {
        "run_accession": args.run,
        "sample": args.sample_description,
        "fastq_files": [str(path) for path in fastqs],
        "read_counts": read_counts,
        "paired_read_count": min(read_counts.values()),
        "junction_supporting_pair_counts": {
            component: len(pairs_by_component[component])
            for component in components
        },
        "junction_supporting_mate_counts": dict(
            Counter(row["component"] for row in evidence_rows)
        ),
        "annotated_transcript_models_matching_each_seed_set": {
            component: sorted(model for model in matching_models[component] if model)
            for component in components
        },
        "interpretation": (
            "Exact junction reads distinguish full-length and predicted "
            "internal DNMT3A transcript structures in developing neocortex."
        ),
        "limit": (
            "An individual library cannot establish developmental regulation; "
            "stage-specific biological replication is required."
        ),
    }
    (RESULTS / f"{args.output_prefix}_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
