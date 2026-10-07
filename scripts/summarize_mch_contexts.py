#!/usr/bin/env python3
"""Classify Bismark CX coverage calls as CG, CA, CC, or CT."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from collections import defaultdict
from pathlib import Path
from typing import TextIO


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXTERNAL_ROOT = Path(
    os.environ.get("DNMT3A_MCH_ROOT", str(ROOT / "external_data/mCH"))
)


class IndexedFasta:
    """Read whole contigs efficiently from an uncompressed FASTA plus .fai."""

    def __init__(self, fasta: Path):
        self.fasta = fasta
        self.index = {}
        with Path(f"{fasta}.fai").open() as handle:
            for line in handle:
                name, length, offset, line_bases, line_width = line.rstrip().split(
                    "\t"
                )[:5]
                self.index[name] = tuple(
                    map(int, (length, offset, line_bases, line_width))
                )
        self.handle = fasta.open("rb")
        self.cached_name: str | None = None
        self.cached_sequence = b""

    def close(self) -> None:
        self.handle.close()

    def contig(self, name: str) -> bytes:
        if name == self.cached_name:
            return self.cached_sequence
        length, offset, line_bases, line_width = self.index[name]
        full_lines, remainder = divmod(length, line_bases)
        byte_count = full_lines * line_width + remainder
        self.handle.seek(offset)
        sequence = self.handle.read(byte_count)
        sequence = sequence.replace(b"\n", b"").replace(b"\r", b"")[:length]
        if len(sequence) != length:
            raise ValueError(f"could not read full FASTA contig {name}")
        self.cached_name = name
        self.cached_sequence = sequence.upper()
        return self.cached_sequence


def oriented_trinucleotide(sequence: bytes, position: int) -> str | None:
    """Return the C-oriented trinucleotide at a 1-based C/G position."""
    index = position - 1
    if not 0 <= index < len(sequence):
        return None
    base = sequence[index : index + 1]
    if base == b"C" and index + 2 < len(sequence):
        return sequence[index : index + 3].decode()
    if base == b"G" and index >= 2:
        complement = bytes.maketrans(b"ACGT", b"TGCA")
        return sequence[index - 2 : index + 1].translate(complement)[::-1].decode()
    return None


def context_labels(trinucleotide: str | None) -> tuple[str, ...]:
    if (
        trinucleotide is None
        or len(trinucleotide) != 3
        or trinucleotide[0] != "C"
        or any(base not in "ACGT" for base in trinucleotide)
    ):
        return ()
    if trinucleotide[1] == "G":
        return ("CG", "all_C")
    dinucleotide = f"C{trinucleotide[1]}"
    broad = "CHG" if trinucleotide[2] == "G" else "CHH"
    return (dinucleotide, broad, "CH", "all_C")


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_gene_bed(
    path: Path,
) -> dict[str, list[tuple[int, int, str, str]]]:
    intervals: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 4:
                raise ValueError("gene BED requires chrom, start, end, and gene")
            strand = fields[5] if len(fields) > 5 else "."
            intervals[fields[0]].append(
                (int(fields[1]), int(fields[2]), fields[3], strand)
            )
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


class GeneSweep:
    """Find gene intervals overlapping monotonically ordered loci."""

    def __init__(self, intervals: dict[str, list[tuple[int, int, str, str]]]):
        self.intervals = intervals
        self.state: dict[str, tuple[int, int, list[tuple[int, int, str, str]]]] = {}

    def overlapping(
        self, chrom: str, position0: int
    ) -> list[tuple[int, int, str, str]]:
        records = self.intervals.get(chrom, [])
        cursor, last, active = self.state.get(chrom, (0, -1, []))
        if position0 < last:
            cursor, active = 0, []
        while cursor < len(records) and records[cursor][0] <= position0:
            active.append(records[cursor])
            cursor += 1
        active = [record for record in active if record[1] > position0]
        self.state[chrom] = (cursor, position0, active)
        return active


def summarize(
    coverage_path: Path,
    fasta_path: Path,
    output_dir: Path,
    bin_size: int = 100_000,
    gene_bed: Path | None = None,
) -> dict[str, object]:
    reference = IndexedFasta(fasta_path)
    totals: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0, 0])
    bins: dict[tuple[str, int, str], list[int]] = defaultdict(lambda: [0, 0, 0])
    gene_intervals = read_gene_bed(gene_bed) if gene_bed else {}
    gene_sweep = GeneSweep(gene_intervals)
    genes: dict[
        tuple[str, int, int, str, str, str], list[int]
    ] = defaultdict(lambda: [0, 0, 0])
    skipped = 0
    rows = 0
    try:
        with open_text(coverage_path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) < 6:
                    raise ValueError(f"invalid coverage row: {line[:100]}")
                chrom = fields[0]
                position = int(fields[1])
                methylated = int(fields[4])
                unmethylated = int(fields[5])
                rows += 1
                if chrom not in reference.index:
                    skipped += 1
                    continue
                labels = context_labels(
                    oriented_trinucleotide(reference.contig(chrom), position)
                )
                if not labels:
                    skipped += 1
                    continue
                contig_class = (
                    "autosome"
                    if chrom.startswith("chr")
                    and chrom[3:].isdigit()
                    and 1 <= int(chrom[3:]) <= 19
                    else "other"
                )
                bin_start = ((position - 1) // bin_size) * bin_size
                for label in labels:
                    for scope in (chrom, contig_class, "genome"):
                        counts = totals[(scope, label)]
                        counts[0] += methylated
                        counts[1] += unmethylated
                        counts[2] += 1
                    if contig_class == "autosome":
                        counts = bins[(chrom, bin_start, label)]
                        counts[0] += methylated
                        counts[1] += unmethylated
                        counts[2] += 1
                    for start, end, gene, strand in gene_sweep.overlapping(
                        chrom, position - 1
                    ):
                        counts = genes[(chrom, start, end, gene, strand, label)]
                        counts[0] += methylated
                        counts[1] += unmethylated
                        counts[2] += 1
    finally:
        reference.close()

    output_dir.mkdir(parents=True, exist_ok=True)
    totals_path = output_dir / "context_summary.tsv"
    with totals_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["scope", "context", "methylated", "unmethylated", "sites", "fraction"]
        )
        for (scope, context), counts in sorted(totals.items()):
            denominator = counts[0] + counts[1]
            fraction = counts[0] / denominator if denominator else float("nan")
            writer.writerow([scope, context, *counts, f"{fraction:.10g}"])

    bins_path = output_dir / "autosomal_bin_context_summary.tsv.gz"
    with gzip.open(bins_path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "chrom",
                "start",
                "end",
                "context",
                "methylated",
                "unmethylated",
                "sites",
                "fraction",
            ]
        )
        for (chrom, start, context), counts in sorted(bins.items()):
            denominator = counts[0] + counts[1]
            fraction = counts[0] / denominator if denominator else float("nan")
            writer.writerow(
                [
                    chrom,
                    start,
                    start + bin_size,
                    context,
                    *counts,
                    f"{fraction:.10g}",
                ]
            )

    gene_path: Path | None = None
    if gene_bed:
        gene_path = output_dir / "gene_body_context_summary.tsv.gz"
        with gzip.open(gene_path, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "chrom",
                    "start",
                    "end",
                    "gene",
                    "strand",
                    "context",
                    "methylated",
                    "unmethylated",
                    "sites",
                    "fraction",
                ]
            )
            for key, counts in sorted(genes.items()):
                denominator = counts[0] + counts[1]
                fraction = counts[0] / denominator if denominator else float("nan")
                writer.writerow([*key, *counts, f"{fraction:.10g}"])

    report = {
        "coverage_file": str(coverage_path),
        "reference": str(fasta_path),
        "input_rows": rows,
        "unclassified_or_unknown_rows": skipped,
        "context_summary": str(totals_path),
        "bin_summary": str(bins_path),
        "gene_bed": str(gene_bed) if gene_bed else None,
        "gene_summary": str(gene_path) if gene_path else None,
    }
    (output_dir / "context_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("coverage", type=Path)
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--gene-bed", type=Path)
    args = parser.parse_args()
    fasta = (
        args.external_root
        / "reference/mm9_controls/mm9_plus_controls.fa"
    )
    print(
        json.dumps(
            summarize(
                args.coverage,
                fasta,
                args.output_dir,
                args.bin_size,
                args.gene_bed,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
