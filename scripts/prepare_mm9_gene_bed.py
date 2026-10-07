#!/usr/bin/env python3
"""Create a longest-transcript protein-coding gene BED from UCSC refGene."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


SOURCE_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz"
)


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def prepare(
    refgene: Path,
    output: Path,
    include_noncoding: bool = False,
) -> int:
    longest: dict[tuple[str, str, str], tuple[int, int, str]] = {}
    with open_text(refgene) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 13:
                continue
            accession = fields[1]
            chrom, strand = fields[2], fields[3]
            start, end = int(fields[4]), int(fields[5])
            gene = fields[12]
            if not include_noncoding and not accession.startswith("NM_"):
                continue
            if not gene or not chrom.startswith("chr") or end <= start:
                continue
            key = (chrom, strand, gene)
            current = longest.get(key)
            if current is None or end - start > current[1] - current[0]:
                longest[key] = (start, end, accession)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for (chrom, strand, gene), (start, end, _) in sorted(
            longest.items(), key=lambda item: (item[0][0], item[1][0], item[1][1])
        ):
            handle.write(f"{chrom}\t{start}\t{end}\t{gene}\t0\t{strand}\n")
    return len(longest)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("refgene", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--include-noncoding", action="store_true")
    args = parser.parse_args()
    count = prepare(args.refgene, args.output, args.include_noncoding)
    print(f"Wrote {count} genes to {args.output}")
    print(f"Source: {SOURCE_URL}")


if __name__ == "__main__":
    main()
