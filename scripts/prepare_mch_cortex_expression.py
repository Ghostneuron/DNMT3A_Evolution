#!/usr/bin/env python3
"""Prepare matched WT cortex expression covariates for gene-body mCA models."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path


DEFAULT_DATA = Path("data/raw/isoform_functional_genomics/GSE164265/cortex_rna")
DEFAULT_MANIFEST = Path("config/cortex_rna_samples.tsv")
EXPERIMENT_COLUMN = {
    "Dnmt3a1_KO": "Dnmt3a1_KO_WT_mean_log2_CPM",
    "Dnmt3a2_KO": "Dnmt3a2_KO_WT_mean_log2_CPM",
    "Dnmt3a1_delta_N": "Dnmt3a1_delta_N_WT_mean_log2_CPM",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensembl_symbol_map(path: Path) -> dict[str, str]:
    mapping = {}
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for identifier in re.findall(r"Ensembl:(ENSMUSG\d+)", row["dbXrefs"]):
                mapping[identifier] = row["Symbol"]
    return mapping


def read_counts(path: Path) -> dict[str, int]:
    counts = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            identifier, count = line.rstrip("\n").split("\t")[:2]
            counts[identifier.split(".")[0]] = int(count)
    return counts


def build(data_dir: Path, manifest_path: Path, output_dir: Path) -> dict[str, object]:
    with manifest_path.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    for sample in manifest:
        path = data_dir / sample["filename"]
        if sha256(path) != sample["sha256"]:
            raise ValueError(f"SHA-256 mismatch for {path}")
    wt_samples = [row for row in manifest if row["genotype"] == "WT"]
    expression: dict[str, dict[str, float]] = defaultdict(dict)
    library_rows = []
    for sample in wt_samples:
        path = data_dir / sample["filename"]
        counts = read_counts(path)
        total = sum(counts.values())
        for gene_id, count in counts.items():
            expression[gene_id][sample["accession"]] = math.log2(
                count * 1_000_000 / total + 1
            )
        library_rows.append(
            {
                "accession": sample["accession"],
                "experiment": sample["experiment"],
                "replicate": sample["replicate"],
                "genes": len(counts),
                "total_counts": total,
                "sha256_valid": True,
            }
        )

    mapping = ensembl_symbol_map(data_dir / "Mus_musculus.gene_info.gz")
    experiments = {
        experiment: [
            row["accession"]
            for row in wt_samples
            if row["experiment"] == experiment
        ]
        for experiment in EXPERIMENT_COLUMN
    }
    rows = []
    for gene_id, sample_values in expression.items():
        symbol = mapping.get(gene_id)
        if symbol is None:
            continue
        row: dict[str, object] = {"ensembl_gene_id": gene_id, "gene": symbol}
        complete = True
        for experiment, column in EXPERIMENT_COLUMN.items():
            values = [sample_values.get(sample) for sample in experiments[experiment]]
            if any(value is None for value in values):
                complete = False
                break
            row[column] = statistics.mean(float(value) for value in values)
        if complete:
            rows.append(row)

    # NCBI occasionally maps retired Ensembl records to the same current symbol.
    # Keep the record with the highest average WT expression deterministically.
    by_symbol: dict[str, dict[str, object]] = {}
    for row in rows:
        score = statistics.mean(float(row[column]) for column in EXPERIMENT_COLUMN.values())
        prior = by_symbol.get(str(row["gene"]))
        if prior is None or score > float(prior["_score"]):
            by_symbol[str(row["gene"])] = {**row, "_score": score}
    rows = [
        {key: value for key, value in row.items() if key != "_score"}
        for _, row in sorted(by_symbol.items())
    ]

    output_dir.mkdir(parents=True, exist_ok=True)
    covariate_path = output_dir / "cortex_WT_expression_covariates.tsv"
    with covariate_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    library_path = output_dir / "cortex_RNA_library_qc.tsv"
    with library_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(library_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(library_rows)

    report = {
        "complete": len(wt_samples) == 8,
        "manifest_libraries_SHA256_validated": len(manifest),
        "WT_libraries": len(wt_samples),
        "mapped_unique_gene_symbols": len(rows),
        "covariate_table": str(covariate_path),
        "library_qc_table": str(library_path),
        "normalization": "log2(counts per million + 1), then mean within matched WT experiment",
        "limitation": "whole cortex expression used as a supporting covariate for NeuN-positive methylomes",
    }
    (output_dir / "cortex_expression_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/mch_gene_body_comparison")
    )
    args = parser.parse_args()
    print(json.dumps(build(args.data_dir, args.manifest, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
