#!/usr/bin/env python3
"""Extract selected DNMT3A records for annotation validation."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/ncbi_dataset/data"
OUT = ROOT / "results/04_selection/mammals/candidate_annotation_validation"
TARGETS = {
    "Carlito_syrichta": {
        "taxname": "Carlito syrichta",
        "gene_id": "103262647",
        "protein": "XP_008058478.1",
        "transcript": "XM_008060287.1",
        "genomic": "NW_007248670.1",
    },
    "Puma_concolor": {
        "taxname": "Puma concolor",
        "gene_id": "112864090",
        "protein": "XP_025782947.1",
        "transcript": "XM_025927162.1",
        "genomic": "NW_020339628.1",
    },
    "Vombatus_ursinus": {
        "taxname": "Vombatus ursinus",
        "gene_id": "114034647",
        "protein": "XP_027706170.1",
        "transcript": "XM_027850369.1",
        "genomic": "NW_020948840.1",
    },
}


def select_records(path: Path, accessions: set[str]):
    selected = []
    for record in SeqIO.parse(path, "fasta"):
        base = record.id.split(":")[0]
        if base in accessions:
            selected.append(record)
    return selected


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reports = {}
    for line in (RAW / "data_report.jsonl").read_text().splitlines():
        report = json.loads(line)
        if report.get("taxname") in {
            values["taxname"] for values in TARGETS.values()
        }:
            reports[report["taxname"]] = report

    rows: list[dict[str, object]] = []
    for safe_id, target in TARGETS.items():
        species_dir = OUT / safe_id
        species_dir.mkdir(parents=True, exist_ok=True)
        specs = {
            "protein": (RAW / "protein.faa", {target["protein"]}, "fasta"),
            "cds": (RAW / "cds.fna", {target["transcript"]}, "fasta"),
            "rna": (RAW / "rna.fna", {target["transcript"]}, "fasta"),
            "gene": (RAW / "gene.fna", {target["genomic"]}, "fasta"),
        }
        lengths = {}
        for label, (source, accessions, fmt) in specs.items():
            records = select_records(source, accessions)
            if len(records) != 1:
                raise RuntimeError(
                    f"Expected one {label} record for {safe_id}, found {len(records)}"
                )
            suffix = "faa" if label == "protein" else "fna"
            output = species_dir / f"{safe_id}.{label}.{suffix}"
            SeqIO.write(records, output, fmt)
            lengths[label] = len(records[0].seq)

        report = reports[target["taxname"]]
        annotation = report["annotations"][0]
        location = annotation["genomicLocations"][0]
        rows.append({
            "safe_id": safe_id,
            "taxname": target["taxname"],
            "gene_id": target["gene_id"],
            "protein_accession": target["protein"],
            "transcript_accession": target["transcript"],
            "genomic_accession": target["genomic"],
            "transcript_type": next(
                item["type"] for item in report["transcriptTypeCounts"]
            ),
            "transcript_count": report["transcriptCount"],
            "protein_count": report["proteinCount"],
            "assembly_accession": annotation["assemblyAccession"],
            "assembly_name": annotation["assemblyName"],
            "annotation_name": annotation["annotationName"],
            "annotation_release_date": annotation["annotationReleaseDate"],
            "scaffold_type": location["sequenceName"],
            "orientation": location["genomicRange"]["orientation"],
            "gene_length_nt": lengths["gene"],
            "rna_length_nt": lengths["rna"],
            "cds_length_nt": lengths["cds"],
            "protein_length_aa": lengths["protein"],
            "rna_minus_cds_nt": lengths["rna"] - lengths["cds"],
        })

    with (OUT / "record_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
