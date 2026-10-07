#!/usr/bin/env python3
"""Summarize matched P21 neuronal Polycomb and DNMT3A bigWigs over mm9 genes."""

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
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_genes(path: Path) -> list[dict[str, object]]:
    genes = {}
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chrom"], int(row["start"]), int(row["end"]), row["gene"], row["strand"])
            genes[key] = {"chrom": key[0], "start": key[1], "end": key[2], "gene": key[3], "strand": key[4]}
    return [genes[key] for key in sorted(genes)]


def mean_signal(bigwig: pyBigWig.pyBigWig, chrom: str, start: int, end: int) -> float:
    chrom_length = bigwig.chroms(chrom)
    if chrom_length is None:
        return math.nan
    start = max(0, min(start, chrom_length))
    end = max(start + 1, min(end, chrom_length))
    value = bigwig.stats(chrom, start, end, type="mean", exact=True)[0]
    return 0.0 if value is None else float(value)


def build(genes_path: Path, manifest_path: Path, data_dir: Path, output_dir: Path, promoter_flank: int) -> dict[str, object]:
    with manifest_path.open(newline="") as handle:
        tracks = list(csv.DictReader(handle, delimiter="\t"))
    for row in tracks:
        path = data_dir / row["filename"]
        if path.stat().st_size != int(row["bytes"]):
            raise ValueError(f"size mismatch for {path}")
        if not row["sha256"] or sha256(path) != row["sha256"]:
            raise ValueError(f"SHA-256 mismatch or missing hash for {path}")
    genes = read_genes(genes_path)
    handles = {row["label"]: pyBigWig.open(str(data_dir / row["filename"])) for row in tracks}
    try:
        results = []
        for gene in genes:
            tss = int(gene["start"]) if gene["strand"] == "+" else int(gene["end"])
            row = {**gene, "length_bp": int(gene["end"]) - int(gene["start"])}
            for label, bigwig in handles.items():
                row[f"{label}_gene_body_mean"] = mean_signal(bigwig, str(gene["chrom"]), int(gene["start"]), int(gene["end"]))
                row[f"{label}_promoter_mean"] = mean_signal(bigwig, str(gene["chrom"]), tss - promoter_flank, tss + promoter_flank)
            results.append(row)
    finally:
        for handle in handles.values():
            handle.close()
    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "p21_polycomb_gene_signals.tsv.gz"
    with gzip.open(table_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    report = {
        "complete": len(tracks) == 6,
        "tracks": len(tracks),
        "genes": len(results),
        "assembly": "mm9",
        "age": "P21",
        "neuronal_marks": "NeuN-positive cortical nuclei",
        "dnmt3a_occupancy": "whole cerebral cortex",
        "promoter_window": f"TSS +/- {promoter_flank} bp",
        "output": str(table_path),
    }
    (output_dir / "p21_polycomb_gene_signals.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", type=Path, default=Path("results/mch_gene_body_comparison/gene_body_mCA_contrasts.tsv.gz"))
    parser.add_argument("--manifest", type=Path, default=Path("config/p21_neuronal_polycomb_tracks.tsv"))
    parser.add_argument("--data-dir", type=Path, default=Path("data/raw/isoform_functional_genomics/GSE164265/p21_polycomb_chromatin"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/mch_gene_body_comparison"))
    parser.add_argument("--promoter-flank", type=int, default=2000)
    args = parser.parse_args()
    print(json.dumps(build(args.genes, args.manifest, args.data_dir, args.output_dir, args.promoter_flank), indent=2))


if __name__ == "__main__":
    main()
