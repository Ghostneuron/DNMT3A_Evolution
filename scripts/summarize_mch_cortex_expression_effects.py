#!/usr/bin/env python3
"""Summarize matched cortex RNA expression effects for DNMT3A perturbations."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

try:
    from scripts.prepare_mch_cortex_expression import (
        DEFAULT_DATA,
        DEFAULT_MANIFEST,
        ensembl_symbol_map,
        read_counts,
        sha256,
    )
except ModuleNotFoundError:
    from prepare_mch_cortex_expression import (
        DEFAULT_DATA,
        DEFAULT_MANIFEST,
        ensembl_symbol_map,
        read_counts,
        sha256,
    )


def build(data_dir: Path, manifest_path: Path, output_dir: Path) -> dict[str, object]:
    with manifest_path.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    normalized: dict[str, dict[str, float]] = defaultdict(dict)
    for sample in manifest:
        path = data_dir / sample["filename"]
        if sha256(path) != sample["sha256"]:
            raise ValueError(f"SHA-256 mismatch for {path}")
        counts = read_counts(path)
        total = sum(counts.values())
        for gene_id, count in counts.items():
            normalized[gene_id][sample["accession"]] = math.log2(
                count * 1_000_000 / total + 1
            )

    mapping = ensembl_symbol_map(data_dir / "Mus_musculus.gene_info.gz")
    rows = []
    for experiment in ("Dnmt3a1_KO", "Dnmt3a2_KO", "Dnmt3a1_delta_N"):
        wt = [
            row["accession"]
            for row in manifest
            if row["experiment"] == experiment and row["genotype"] == "WT"
        ]
        treatment = [
            row["accession"]
            for row in manifest
            if row["experiment"] == experiment and row["genotype"] == experiment
        ]
        for gene_id, values in normalized.items():
            symbol = mapping.get(gene_id)
            if symbol is None or any(sample not in values for sample in wt + treatment):
                continue
            wt_values = [values[sample] for sample in wt]
            treatment_values = [values[sample] for sample in treatment]
            rows.append(
                {
                    "ensembl_gene_id": gene_id,
                    "gene": symbol,
                    "contrast": f"{experiment}-WT",
                    "wt_n": len(wt),
                    "treatment_n": len(treatment),
                    "wt_mean_log2_CPM": statistics.mean(wt_values),
                    "treatment_mean_log2_CPM": statistics.mean(treatment_values),
                    "expression_difference": statistics.mean(treatment_values)
                    - statistics.mean(wt_values),
                }
            )

    unique = {}
    for row in rows:
        key = (row["contrast"], row["gene"])
        prior = unique.get(key)
        if prior is None or float(row["wt_mean_log2_CPM"]) > float(
            prior["wt_mean_log2_CPM"]
        ):
            unique[key] = row
    rows = [unique[key] for key in sorted(unique)]

    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "cortex_expression_contrasts.tsv.gz"
    import gzip

    with gzip.open(table_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    counts_by_contrast = {
        contrast: sum(row["contrast"] == contrast for row in rows)
        for contrast in sorted({str(row["contrast"]) for row in rows})
    }
    report = {
        "complete": len(manifest) == 16,
        "SHA256_validated_libraries": len(manifest),
        "genes_by_contrast": counts_by_contrast,
        "output": str(table_path),
        "effect": "mean treatment log2(CPM + 1) minus matched WT mean",
        "warning": "descriptive whole-cortex expression effect; not a formal differential-expression test",
    }
    (output_dir / "cortex_expression_contrasts.json").write_text(
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
