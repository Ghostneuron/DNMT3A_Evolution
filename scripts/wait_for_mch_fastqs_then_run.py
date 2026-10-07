#!/usr/bin/env python3
"""Wait for checksum-valid FASTQs, then launch one full EM-seq sample."""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
import time
from pathlib import Path

try:
    from scripts.run_mch_bismark import DEFAULT_EXTERNAL_ROOT, MANIFEST
except ModuleNotFoundError:
    from run_mch_bismark import DEFAULT_EXTERNAL_ROOT, MANIFEST


ROOT = Path(__file__).resolve().parents[1]
PYTHON = Path("/opt/anaconda3/bin/python")
PIPELINE = ROOT / "scripts/run_mch_bismark.py"


def md5(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest_files() -> dict[str, tuple[int, str]]:
    with MANIFEST.open(newline="") as handle:
        return {
            Path(row["fastq_url"]).name: (int(row["bytes"]), row["md5"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--wait-for-file", action="append", required=True)
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--poll-seconds", type=int, default=120)
    args = parser.parse_args()

    expected = manifest_files()
    raw = args.external_root / "raw/GSE164265/fastq"
    required = [raw / name for name in args.wait_for_file]
    for path in required:
        if path.name not in expected:
            raise RuntimeError(f"file is absent from manifest: {path.name}")

    while True:
        ready = True
        for path in required:
            size, _ = expected[path.name]
            if (
                not path.is_file()
                or path.stat().st_size != size
                or Path(f"{path}.aria2").exists()
            ):
                ready = False
                break
        if ready:
            break
        time.sleep(args.poll_seconds)

    for path in required:
        _, checksum = expected[path.name]
        observed = md5(path)
        if observed != checksum:
            raise RuntimeError(
                f"checksum mismatch for {path.name}: {observed} != {checksum}"
            )
        print(f"validated {path.name}", flush=True)

    command = [
        str(PYTHON),
        str(PIPELINE),
        "all",
        "--sample",
        args.sample,
        "--threads",
        str(args.threads),
        "--ignore",
        "10",
        "--ignore-r2",
        "5",
        "--ignore-3prime",
        "10",
        "--ignore-3prime-r2",
        "0",
    ]
    subprocess.run(command, check=True, cwd=ROOT)


if __name__ == "__main__":
    main()
