#!/usr/bin/env python3
"""Annotate top GSE164265 isoform-differential windows with mm9 RefSeq genes."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFGENE = ROOT / "data/raw/mm9_annotation/refGene.txt.gz"
INPUT = (
    ROOT
    / "results/07_isoform_function/GSE164265_neuron_nuclei_methylation"
    / "top_isoform_differential_windows.tsv"
)
OUTPUT = (
    ROOT
    / "results/07_isoform_function/GSE164265_neuron_nuclei_methylation"
    / "top_isoform_differential_windows_annotated.tsv"
)


def interval_distance(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    if end_a <= start_b:
        return start_b - end_a
    if end_b <= start_a:
        return start_a - end_b
    return 0


def read_genes(path: Path = REFGENE) -> dict[str, list[tuple[int, int, str]]]:
    collapsed: dict[tuple[str, str], list[int]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 13:
                continue
            chromosome = fields[2]
            start, end = int(fields[4]), int(fields[5])
            symbol = fields[12]
            key = (chromosome, symbol)
            if key not in collapsed:
                collapsed[key] = [start, end]
            else:
                collapsed[key][0] = min(collapsed[key][0], start)
                collapsed[key][1] = max(collapsed[key][1], end)
    genes: defaultdict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for (chromosome, symbol), (start, end) in collapsed.items():
        genes[chromosome].append((start, end, symbol))
    for chromosome in genes:
        genes[chromosome].sort()
    return dict(genes)


def annotate_window(
    chromosome: str,
    start: int,
    end: int,
    genes: dict[str, list[tuple[int, int, str]]],
) -> tuple[str, str, int | str]:
    chromosome_genes = genes.get(chromosome, [])
    overlaps = sorted({
        symbol
        for gene_start, gene_end, symbol in chromosome_genes
        if interval_distance(start, end, gene_start, gene_end) == 0
    })
    distances = [
        (interval_distance(start, end, gene_start, gene_end), symbol)
        for gene_start, gene_end, symbol in chromosome_genes
    ]
    if not distances:
        return ";".join(overlaps), "", ""
    minimum = min(distance for distance, _ in distances)
    nearest = sorted({symbol for distance, symbol in distances if distance == minimum})
    return ";".join(overlaps), ";".join(nearest), minimum


def main() -> None:
    genes = read_genes()
    with INPUT.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    output = []
    for row in rows:
        overlaps, nearest, distance = annotate_window(
            row["chromosome"], int(row["start"]), int(row["end"]), genes
        )
        output.append({
            **row,
            "overlapping_RefSeq_genes": overlaps,
            "nearest_RefSeq_genes": nearest,
            "nearest_gene_distance_bp": distance,
        })
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(output[0])
        )
        writer.writeheader()
        writer.writerows(output)
    print(f"annotated {len(output)} windows")


if __name__ == "__main__":
    main()
