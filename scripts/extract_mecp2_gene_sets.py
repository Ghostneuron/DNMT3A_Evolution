#!/usr/bin/env python3
"""Extract independent MeCP2 neuronal gene classes from Moore et al. 2025."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import openpyxl


EXPECTED_SHA256 = "afec580bdf173ada6fb940b8d8f5517905cfe21832957a223165238ddc120461"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build(input_path: Path, output_dir: Path) -> dict[str, object]:
    if sha256(input_path) != EXPECTED_SHA256:
        raise ValueError("MeCP2 supplementary workbook SHA-256 mismatch")
    workbook = openpyxl.load_workbook(input_path, read_only=True, data_only=True)
    sheet = workbook["TableS2_subclass_core_label"]
    iterator = sheet.iter_rows(values_only=True)
    header = [str(value) for value in next(iterator)]
    rows = []
    for values in iterator:
        row = dict(zip(header, values))
        subclasses = [str(row[name]) for name in ("L4", "L5", "PV", "SST")]
        rows.append(
            {
                "gene": row["gene"],
                "L4": row["L4"],
                "L5": row["L5"],
                "PV": row["PV"],
                "SST": row["SST"],
                "Core": row["Core"],
                "MR_subclass_count": sum(value == "MR" for value in subclasses),
                "MA_subclass_count": sum(value == "MA" for value in subclasses),
                "all_subclasses_unchanged": all(
                    value == "Unchanged" for value in subclasses
                ),
            }
        )
    workbook.close()

    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "moore_2025_mecp2_gene_classes.tsv"
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "complete": True,
        "source_sheet": "TableS2_subclass_core_label",
        "genes": len(rows),
        "core_MR": sum(row["Core"] == "MR" for row in rows),
        "core_MA": sum(row["Core"] == "MA" for row in rows),
        "core_unchanged": sum(row["Core"] == "Unchanged" for row in rows),
        "output": str(table_path),
    }
    (output_dir / "moore_2025_mecp2_gene_classes.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/raw/mecp2_gene_sets/Moore_2025_Supplemental_Tables.xlsx"),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/mch_gene_body_comparison")
    )
    args = parser.parse_args()
    print(json.dumps(build(args.input, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
