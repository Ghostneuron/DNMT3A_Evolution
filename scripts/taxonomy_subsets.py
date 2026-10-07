#!/usr/bin/env python3
"""Build reproducible NCBI-taxonomy metadata and DNMT3A alignment subsets."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPRESENTATIVES = ROOT / "results/01_curated/representatives.tsv"
DATA_REPORT = ROOT / "data/raw/ncbi_dataset/data/data_report.jsonl"
TAXONOMY_IDS = ROOT / "data/provenance/ncbi_taxonomy_ids.txt"
TAXONOMY_XML = ROOT / "data/provenance/ncbi_taxonomy_efetch_2026-07-22.xml"
TAXON_METADATA = ROOT / "results/01_curated/taxon_metadata.tsv"
SEQUENCE_QC = ROOT / "results/02_alignment/sequence_codon_qc.tsv"
ALIGNMENT = ROOT / "results/02_alignment/DNMT3A_MACSE_clean_HyPhy.fasta"
SUBSET_DIR = ROOT / "results/02_alignment/subsets"
CONSTRAINT_DIR = ROOT / "results/03_tree/constraints"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(chunks)))
                identifier, chunks = line[1:].split()[0], []
            else:
                if identifier is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if identifier is not None:
        records.append((identifier, "".join(chunks)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start:start + width] + "\n")


def selected_species_taxids() -> list[dict[str, str]]:
    if not REPRESENTATIVES.exists() or not DATA_REPORT.exists():
        raise FileNotFoundError("Run DNMT3A curation and retain the NCBI data report first")
    representatives = read_tsv(REPRESENTATIVES)
    reports = [json.loads(line) for line in DATA_REPORT.read_text().splitlines() if line.strip()]
    report_by_name = {row["taxname"]: row for row in reports}
    rows: list[dict[str, str]] = []
    for representative in representatives:
        species = representative["species"]
        if species not in report_by_name:
            raise ValueError(f"No NCBI gene-report taxon match for {species}")
        rows.append({
            "safe_id": representative["safe_id"],
            "species": species,
            "taxid": str(report_by_name[species]["taxId"]),
        })
    if len({row["taxid"] for row in rows}) != len(rows):
        raise ValueError("Selected representatives do not map one-to-one to NCBI taxon IDs")
    return rows


def write_ids() -> None:
    rows = selected_species_taxids()
    TAXONOMY_IDS.parent.mkdir(parents=True, exist_ok=True)
    TAXONOMY_IDS.write_text(",".join(sorted((row["taxid"] for row in rows), key=int)) + "\n")
    print(f"Wrote {len(rows)} taxon IDs to {TAXONOMY_IDS}")


def taxonomy_records() -> dict[str, dict[str, object]]:
    if not TAXONOMY_XML.exists():
        raise FileNotFoundError(
            f"Missing {TAXONOMY_XML}; fetch it with the documented NCBI EFetch command"
        )
    root = ET.parse(TAXONOMY_XML).getroot()
    records: dict[str, dict[str, object]] = {}
    for taxon in root.findall("./Taxon"):
        taxid = taxon.findtext("TaxId", "")
        lineage_nodes = taxon.findall("./LineageEx/Taxon")
        nodes = [
            {
                "taxid": node.findtext("TaxId", ""),
                "name": node.findtext("ScientificName", ""),
                "rank": node.findtext("Rank", ""),
            }
            for node in lineage_nodes
        ]
        nodes.append({
            "taxid": taxid,
            "name": taxon.findtext("ScientificName", ""),
            "rank": taxon.findtext("Rank", ""),
        })
        records[taxid] = {
            "scientific_name": taxon.findtext("ScientificName", ""),
            "nodes": nodes,
        }
    return records


def classify(names: set[str]) -> str:
    for clade in ("Mammalia", "Aves", "Crocodylia", "Testudines", "Squamata", "Amphibia", "Dipnoi"):
        if clade in names:
            return clade
    return "other_sarcopterygian"


def build() -> None:
    selected = selected_species_taxids()
    taxonomy = taxonomy_records()
    missing = sorted({row["taxid"] for row in selected} - set(taxonomy), key=int)
    if missing:
        raise ValueError(f"NCBI taxonomy response is missing {len(missing)} selected IDs: {missing[:10]}")

    metadata: list[dict[str, object]] = []
    for row in selected:
        record = taxonomy[row["taxid"]]
        nodes = record["nodes"]
        assert isinstance(nodes, list)
        ranks = {
            str(node["rank"]): str(node["name"])
            for node in nodes
            if node["rank"] and node["rank"] != "no rank"
        }
        names = {str(node["name"]) for node in nodes}
        major_clade = classify(names)
        metadata.append({
            **row,
            "ncbi_scientific_name": record["scientific_name"],
            "major_clade": major_clade,
            "is_mammal": str("Mammalia" in names).lower(),
            "class": ranks.get("class", ""),
            "order": ranks.get("order", ""),
            "family": ranks.get("family", ""),
            "genus": ranks.get("genus", ""),
            "lineage": ";".join(str(node["name"]) for node in nodes),
        })
    fields = [
        "safe_id", "species", "taxid", "ncbi_scientific_name", "major_clade",
        "is_mammal", "class", "order", "family", "genus", "lineage",
    ]
    write_tsv(TAXON_METADATA, metadata, fields)

    qc = {row["safe_id"]: row for row in read_tsv(SEQUENCE_QC)}
    alignment = read_fasta(ALIGNMENT)
    metadata_by_id = {str(row["safe_id"]): row for row in metadata}
    alignment_ids = {identifier for identifier, _ in alignment}
    if alignment_ids != set(metadata_by_id) or alignment_ids != set(qc):
        raise ValueError("Taxonomy, sequence-QC, and alignment IDs are not identical")

    sequence_groups: dict[str, list[str]] = defaultdict(list)
    for identifier, sequence in alignment:
        sequence_groups[sequence].append(identifier)
    identical_rows: list[dict[str, object]] = []
    for group_number, members in enumerate(
        sorted((members for members in sequence_groups.values() if len(members) > 1), key=lambda x: x[0]),
        start=1,
    ):
        identical_rows.append({
            "group": f"identical_{group_number}",
            "sequences": len(members),
            "safe_ids": ",".join(members),
            "species": ";".join(str(metadata_by_id[member]["species"]) for member in members),
        })
    write_tsv(
        SUBSET_DIR / "identical_sequence_groups.tsv",
        identical_rows,
        ["group", "sequences", "safe_ids", "species"],
    )

    high_coverage = {
        identifier
        for identifier, row in qc.items()
        if float(row["retained_fraction_after_column_filter"]) >= 0.85
    }
    mammals = {
        identifier for identifier, row in metadata_by_id.items()
        if row["is_mammal"] == "true"
    }
    subset_ids = {
        "all_475": alignment_ids,
        "all_high_coverage": high_coverage,
        "mammals": mammals,
        "mammals_high_coverage": mammals & high_coverage,
    }
    manifest: list[dict[str, object]] = []
    for name, identifiers in subset_ids.items():
        output = SUBSET_DIR / f"DNMT3A_{name}.fasta"
        records = [(identifier, sequence) for identifier, sequence in alignment if identifier in identifiers]
        write_fasta(output, records)
        manifest.append({
            "subset": name,
            "sequences": len(records),
            "codon_columns": len(records[0][1]) // 3 if records else 0,
            "minimum_sequence_coverage": "0.85" if "high_coverage" in name else "none",
            "path": str(output.relative_to(ROOT)),
        })
    write_tsv(
        SUBSET_DIR / "subset_manifest.tsv",
        manifest,
        ["subset", "sequences", "codon_columns", "minimum_sequence_coverage", "path"],
    )

    mammal_orders: dict[str, list[str]] = defaultdict(list)
    for identifier in sorted(mammals):
        order = str(metadata_by_id[identifier]["order"])
        if not order:
            raise ValueError(f"Mammalian representative lacks NCBI order: {identifier}")
        mammal_orders[order].append(identifier)
    order_clades = [
        f"({','.join(sorted(members))}){order}" if len(members) > 1 else members[0]
        for order, members in sorted(mammal_orders.items())
    ]
    CONSTRAINT_DIR.mkdir(parents=True, exist_ok=True)
    (CONSTRAINT_DIR / "NCBI_mammal_orders_constraint.nwk").write_text(
        f"({','.join(order_clades)})Mammalia;\n"
    )
    counts: dict[str, int] = {}
    for row in metadata:
        clade = str(row["major_clade"])
        counts[clade] = counts.get(clade, 0) + 1
    print(json.dumps({
        "taxonomy_records": len(metadata),
        "major_clades": dict(sorted(counts.items())),
        "subsets": {row["subset"]: row["sequences"] for row in manifest},
        "identical_sequence_groups": len(identical_rows),
        "mammalian_orders_in_constraint": len(mammal_orders),
    }, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["ids", "build"])
    args = parser.parse_args()
    if args.stage == "ids":
        write_ids()
    else:
        build()


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}")
