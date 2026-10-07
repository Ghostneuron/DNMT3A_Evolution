#!/usr/bin/env python3
"""Combine source-specific phenotype tables into the validated project table."""

from __future__ import annotations

import csv
from pathlib import Path

try:
    from .import_brainzoo_mch import FIELDS, ROOT
except ImportError:  # Direct script execution.
    from import_brainzoo_mch import FIELDS, ROOT


SOURCE_DIR = ROOT / "data/traits/sources"
OUTPUT = ROOT / "data/traits/brain_phenotypes.tsv"


def main() -> None:
    rows: list[dict[str, str]] = []
    paths = sorted(SOURCE_DIR.glob("*.tsv"))
    if not paths:
        raise SystemExit(f"No phenotype source tables found in {SOURCE_DIR}")
    for path in paths:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != FIELDS:
                raise SystemExit(f"{path}: phenotype schema mismatch")
            rows.extend(reader)
    rows.sort(key=lambda row: (
        row["species"], row["phenotype"], row["developmental_stage"],
        row["brain_region"], row["sample_id"],
    ))
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} records from {len(paths)} source tables to {OUTPUT}")


if __name__ == "__main__":
    main()
