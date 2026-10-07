#!/usr/bin/env python3
"""Restore and checksum-validate the compact-analysis source bigWigs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(rows: list[dict[str, str]], output_dir: Path) -> list[dict[str, object]]:
    results = []
    for row in rows:
        path = output_dir / row["filename"]
        observed_bytes = path.stat().st_size if path.exists() else 0
        observed_hash = sha256(path) if observed_bytes == int(row["bytes"]) else ""
        results.append(
            {
                "accession": row["accession"],
                "path": str(path),
                "expected_bytes": int(row["bytes"]),
                "observed_bytes": observed_bytes,
                "size_matches": observed_bytes == int(row["bytes"]),
                "sha256_matches": observed_hash == row["sha256"],
            }
        )
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("config/cortex_chromatin_tracks.tsv"))
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("data/raw/isoform_functional_genomics/GSE164265/cortex_chromatin"),
    )
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if not args.validate_only:
        for row in rows:
            subprocess.run(
                ["curl", "--fail", "--location", "--continue-at", "-", "--output", str(args.output_dir / row["filename"]), row["url"]],
                check=True,
            )
    results = validate(rows, args.output_dir)
    failed = [row for row in results if not row["sha256_matches"]]
    print(json.dumps({"files": len(results), "validated": len(results) - len(failed), "details": results}, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
