#!/usr/bin/env python3
"""Download and verify the eight GSE164265 neuronal EM-seq libraries."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config/mch_samples.tsv"
DEFAULT_EXTERNAL_ROOT = Path(
    os.environ.get("DNMT3A_MCH_ROOT", str(ROOT / "external_data/mCH"))
)


def md5(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest(path: Path = MANIFEST) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def aria2_input_lines(rows: list[dict[str, str]]) -> list[str]:
    """Build one checksum-enforced aria2 input block per FASTQ."""
    lines = []
    for row in rows:
        lines.extend([
            row["fastq_url"],
            f"  out={Path(row['fastq_url']).name}",
            f"  checksum=md5={row['md5']}",
        ])
    return lines


def validate_completed(
    rows: list[dict[str, str]],
    output_dir: Path,
) -> list[dict[str, str | int | bool]]:
    results = []
    for row in rows:
        filename = Path(row["fastq_url"]).name
        path = output_dir / filename
        expected_size = int(row["bytes"])
        size = path.stat().st_size if path.exists() else 0
        size_matches = size == expected_size
        checksum = md5(path) if size_matches else ""
        results.append({
            "run_accession": row["run_accession"],
            "read_pair": int(row["read_pair"]),
            "local_file": str(path),
            "expected_bytes": expected_size,
            "observed_bytes": size,
            "size_matches": size_matches,
            "expected_md5": row["md5"],
            "observed_md5": checksum,
            "md5_matches": checksum == row["md5"] if checksum else False,
        })
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT
    )
    parser.add_argument("--connections", type=int, default=4)
    parser.add_argument(
        "--connections-per-file",
        type=int,
        default=4,
        help="Number of aria2 segments/connections per FASTQ (use 1 for repair)",
    )
    parser.add_argument(
        "--sample",
        action="append",
        help="Restrict to one SRR run accession; repeat to select several",
    )
    parser.add_argument(
        "--file",
        action="append",
        dest="filename",
        help="Restrict to one FASTQ filename; repeat to select several",
    )
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if args.connections < 1:
        parser.error("--connections must be at least 1")
    if args.connections_per_file < 1:
        parser.error("--connections-per-file must be at least 1")

    rows = read_manifest()
    if args.sample:
        selected = set(args.sample)
        rows = [row for row in rows if row["run_accession"] in selected]
        missing = selected - {row["run_accession"] for row in rows}
        if missing:
            parser.error(f"unknown run accession(s): {', '.join(sorted(missing))}")
    if args.filename:
        selected_files = set(args.filename)
        rows = [
            row for row in rows
            if Path(row["fastq_url"]).name in selected_files
        ]
        missing_files = selected_files - {
            Path(row["fastq_url"]).name for row in rows
        }
        if missing_files:
            parser.error(
                f"unknown FASTQ filename(s): {', '.join(sorted(missing_files))}"
            )
    output_dir = args.external_root / "raw/GSE164265/fastq"
    output_dir.mkdir(parents=True, exist_ok=True)
    if not args.validate_only:
        input_path = args.external_root / "logs/aria2_fastq_input.txt"
        input_path.parent.mkdir(parents=True, exist_ok=True)
        input_lines = aria2_input_lines(rows)
        input_path.write_text("\n".join(input_lines) + "\n")
        command = [
            "/opt/homebrew/bin/aria2c",
            "--continue=true",
            "--auto-file-renaming=false",
            "--allow-overwrite=false",
            f"--max-concurrent-downloads={args.connections}",
            f"--split={args.connections_per_file}",
            f"--max-connection-per-server={args.connections_per_file}",
            "--min-split-size=64M",
            "--file-allocation=none",
            "--show-console-readout=false",
            "--summary-interval=60",
            f"--dir={output_dir}",
            f"--input-file={input_path}",
        ]
        subprocess.run(command, check=True)

    validation = validate_completed(rows, output_dir)
    validation_name = (
        "fastq_validation.json"
        if not args.sample and not args.filename
        else "fastq_validation."
        + ".".join(sorted(args.sample or []) + sorted(args.filename or []))
        + ".json"
    )
    validation_path = args.external_root / "raw/GSE164265" / validation_name
    validation_path.write_text(json.dumps(validation, indent=2) + "\n")
    failed = [row for row in validation if not row["md5_matches"]]
    print(json.dumps({
        "files": len(validation),
        "validated": len(validation) - len(failed),
        "failed_or_incomplete": len(failed),
        "validation_file": str(validation_path),
    }, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
