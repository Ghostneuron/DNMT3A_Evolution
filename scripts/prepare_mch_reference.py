#!/usr/bin/env python3
"""Create an mm9 plus NEB EM-seq control composite reference."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXTERNAL_ROOT = Path(
    os.environ.get("DNMT3A_MCH_ROOT", str(ROOT / "external_data/mCH"))
)
CONTROL_SOURCE = (
    "https://raw.githubusercontent.com/nebiolabs/EM-seq/"
    "3f51e84f5e841ffd560513b1c99943bb8148bf0b/"
    "assets/methylation_controls.fa"
)


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_headers(path: Path) -> list[str]:
    headers = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def build_composite(mouse_fasta: Path, control_fasta: Path, output: Path) -> None:
    controls = fasta_headers(control_fasta)
    normalized = [header.lower() for header in controls]
    has_lambda = any(header.endswith("lambda") for header in normalized)
    has_puc19 = any(header.endswith("puc19c") for header in normalized)
    if not has_lambda or not has_puc19:
        raise ValueError(
            "control FASTA must contain NEB lambda and pUC19c sequences"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as destination:
        with mouse_fasta.open("rb") as source:
            shutil.copyfileobj(source, destination, 8 * 1024 * 1024)
        destination.write(b"\n")
        with control_fasta.open("rb") as source:
            shutil.copyfileobj(source, destination, 1024 * 1024)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    parser.add_argument(
        "--control-fasta",
        type=Path,
        help="Downloaded NEB methylation_controls.fa",
    )
    args = parser.parse_args()

    mouse_fasta = args.external_root / "reference/mm9/mm9.fa"
    control_fasta = args.control_fasta or (
        args.external_root / "reference/controls/methylation_controls.fa"
    )
    output = (
        args.external_root / "reference/mm9_controls/mm9_plus_controls.fa"
    )
    build_composite(mouse_fasta, control_fasta, output)
    provenance = {
        "mouse_fasta": str(mouse_fasta),
        "mouse_sha256": sha256(mouse_fasta),
        "control_fasta": str(control_fasta),
        "control_sha256": sha256(control_fasta),
        "control_source": CONTROL_SOURCE,
        "control_sequences": fasta_headers(control_fasta),
        "composite_fasta": str(output),
        "composite_sha256": sha256(output),
    }
    provenance_path = output.parent / "reference_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
