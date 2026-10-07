#!/usr/bin/env python3
"""Safely continue one EM-seq sample after an externally running alignment.

The script waits for the alignment wrapper PID, validates each BAM before
removing its superseded predecessor, and then runs deduplication and methylation
extraction through the project's normal pipeline wrapper.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PIPELINE = ROOT / "scripts/run_mch_bismark.py"
PYTHON = Path("/opt/anaconda3/bin/python")
SAMTOOLS = Path("/opt/homebrew/bin/samtools")
DEFAULT_EXTERNAL_ROOT = Path(
    os.environ.get("DNMT3A_MCH_ROOT", str(ROOT / "external_data/mCH"))
)


def process_exists(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def validate_bam(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"missing or empty BAM: {path}")
    subprocess.run([str(SAMTOOLS), "quickcheck", "-v", str(path)], check=True)


def remove_files(directory: Path) -> int:
    reclaimed = 0
    if not directory.is_dir():
        return reclaimed
    for path in directory.iterdir():
        if path.is_file():
            reclaimed += path.stat().st_size
            path.unlink()
    try:
        directory.rmdir()
    except OSError:
        pass
    return reclaimed


def run_stage(
    stage: str,
    sample: str,
    external_root: Path,
    threads: int,
    extraction_trimming: tuple[int, int, int, int] | None = None,
) -> None:
    command = [
        str(PYTHON),
        str(PIPELINE),
        stage,
        "--external-root",
        str(external_root),
        "--sample",
        sample,
        "--threads",
        str(threads),
    ]
    if extraction_trimming is not None:
        command.extend(
            [
                "--ignore",
                str(extraction_trimming[0]),
                "--ignore-r2",
                str(extraction_trimming[1]),
                "--ignore-3prime",
                str(extraction_trimming[2]),
                "--ignore-3prime-r2",
                str(extraction_trimming[3]),
            ]
        )
    subprocess.run(command, check=True, cwd=ROOT)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--alignment-pid", required=True, type=int)
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()

    sample_root = args.external_root / "results/bismark" / args.sample
    aligned_bam = sample_root / "alignment" / f"{args.sample}_pe.bam"
    alignment_report = (
        sample_root / "alignment" / f"{args.sample}_PE_report.txt"
    )
    deduplicated_bam = (
        sample_root / "deduplicated" / f"{args.sample}_pe.deduplicated.bam"
    )

    print(f"Waiting for alignment wrapper PID {args.alignment_pid}", flush=True)
    while process_exists(args.alignment_pid):
        time.sleep(args.poll_seconds)

    if not alignment_report.is_file() or alignment_report.stat().st_size == 0:
        raise RuntimeError(f"alignment report was not finalized: {alignment_report}")
    validate_bam(aligned_bam)
    reclaimed = remove_files(sample_root / "alignment" / "tmp")
    print(f"Validated alignment; reclaimed {reclaimed / 2**30:.1f} GiB", flush=True)

    run_stage("dedup", args.sample, args.external_root, args.threads)
    validate_bam(deduplicated_bam)
    aligned_size = aligned_bam.stat().st_size
    aligned_bam.unlink()
    print(
        f"Validated deduplicated BAM; reclaimed {aligned_size / 2**30:.1f} GiB",
        flush=True,
    )

    run_stage(
        "extract",
        args.sample,
        args.external_root,
        args.threads,
        extraction_trimming=(10, 5, 10, 0),
    )
    print("Methylation extraction completed", flush=True)


if __name__ == "__main__":
    main()
