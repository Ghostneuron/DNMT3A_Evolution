#!/usr/bin/env python3
"""Summarize source-study cortex bigWig tracks over analyzed mm9 genes."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import pyBigWig


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest(path: Path, data_dir: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        source = data_dir / row["filename"]
        if source.stat().st_size != int(row["bytes"]):
            raise ValueError(f"size mismatch for {source}")
        if row["sha256"] and sha256(source) != row["sha256"]:
            raise ValueError(f"SHA-256 mismatch for {source}")
    return rows


def read_genes(path: Path) -> list[dict[str, object]]:
    genes = {}
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chrom"], int(row["start"]), int(row["end"]), row["gene"], row["strand"])
            genes[key] = {
                "chrom": key[0], "start": key[1], "end": key[2],
                "gene": key[3], "strand": key[4],
            }
    return [genes[key] for key in sorted(genes)]


def mean_signal(bigwig: pyBigWig.pyBigWig, chrom: str, start: int, end: int) -> float:
    chrom_length = bigwig.chroms(chrom)
    if chrom_length is None:
        return math.nan
    start = max(0, min(start, chrom_length))
    end = max(start + 1, min(end, chrom_length))
    value = bigwig.stats(chrom, start, end, type="mean", exact=True)[0]
    return 0.0 if value is None else float(value)


def build(
    genes_path: Path,
    manifest_path: Path,
    data_dir: Path,
    output_dir: Path,
    promoter_flank: int,
) -> dict[str, object]:
    tracks = read_manifest(manifest_path, data_dir)
    genes = read_genes(genes_path)
    handles = {
        row["label"]: pyBigWig.open(str(data_dir / row["filename"])) for row in tracks
    }
    try:
        rows = []
        for gene in genes:
            tss = int(gene["start"]) if gene["strand"] == "+" else int(gene["end"])
            row = {**gene, "length_bp": int(gene["end"]) - int(gene["start"])}
            for label, bigwig in handles.items():
                row[f"{label}_gene_body_mean"] = mean_signal(
                    bigwig, str(gene["chrom"]), int(gene["start"]), int(gene["end"])
                )
                row[f"{label}_promoter_mean"] = mean_signal(
                    bigwig, str(gene["chrom"]), tss - promoter_flank, tss + promoter_flank
                )
            rows.append(row)
    finally:
        for handle in handles.values():
            handle.close()

    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "cortex_chromatin_gene_signals.tsv.gz"
    with gzip.open(table_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "complete": len(tracks) == 4,
        "tracks": len(tracks),
        "genes": len(rows),
        "assembly": "mm9",
        "tissue_age": "P18 cerebral cortex",
        "promoter_window": f"TSS +/- {promoter_flank} bp",
        "output": str(table_path),
    }
    (output_dir / "cortex_chromatin_gene_signals.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genes", type=Path,
        default=Path("results/mch_gene_body_comparison/gene_body_mCA_contrasts.tsv.gz"),
    )
    parser.add_argument(
        "--manifest", type=Path, default=Path("config/cortex_chromatin_tracks.tsv")
    )
    parser.add_argument(
        "--data-dir", type=Path,
        default=Path("data/raw/isoform_functional_genomics/GSE164265/cortex_chromatin"),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/mch_gene_body_comparison")
    )
    parser.add_argument("--promoter-flank", type=int, default=2000)
    args = parser.parse_args()
    print(json.dumps(build(args.genes, args.manifest, args.data_dir, args.output_dir, args.promoter_flank), indent=2))


if __name__ == "__main__":
    main()
