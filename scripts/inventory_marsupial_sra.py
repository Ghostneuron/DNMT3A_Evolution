#!/usr/bin/env python3
"""Build a reproducible inventory of marsupial SRA datasets relevant to DNMT3A.

The inventory deliberately separates evidence for an exact capped TSS from
evidence for a transcript first exon. Ordinary RNA-seq and cDNA long reads
cannot, by themselves, establish the exact transcription-initiation base.
"""

from __future__ import annotations

import csv
import json
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/marsupial_sra"
RESULTS = ROOT / "results/06_isoform_evolution"
INPUTS = {
    "fat-tailed dunnart": RAW / "dunnart_sra_records.xml",
    "tammar wallaby": RAW / "tammar_sra_records.xml",
}

BRAIN_TERMS = (
    "brain", "cerebell", "cortex", "cortical", "neocortex", "neural",
    "forebrain", "hindbrain", "hippocamp", "striat", "thalam",
)
DEVELOPMENT_TERMS = (
    "embry", "fetal", "foetal", "pouch young", "neonat", "postnatal",
    "development", "p12", "p20",
)
CAP_STRATEGIES = ("cage", "rampage", "5' race", "5-prime race", "5 prime race")
LONG_READ_PLATFORMS = ("OXFORD_NANOPORE", "PACBIO_SMRT")
TRANSCRIPTOMIC_STRATEGIES = (
    "RNA-Seq", "EST", "FL-cDNA", "OTHER",
)


def text(node: ET.Element | None, path: str, default: str = "") -> str:
    if node is None:
        return default
    value = node.findtext(path)
    return value.strip() if value else default


def external_id(node: ET.Element | None, namespace: str) -> str:
    if node is None:
        return ""
    for item in node.findall(".//EXTERNAL_ID"):
        if item.get("namespace") == namespace:
            return (item.text or "").strip()
    return ""


def attributes(node: ET.Element | None, container: str, item: str) -> dict[str, str]:
    values: dict[str, str] = {}
    if node is None:
        return values
    for attribute in node.findall(f"./{container}/{item}"):
        tag = text(attribute, "TAG").lower()
        value = text(attribute, "VALUE")
        if tag and value:
            values[tag] = value
    return values


def classify_evidence(
    strategy: str,
    selection: str,
    platform: str,
    searchable_text: str,
) -> tuple[str, int, str]:
    lower = searchable_text.lower()
    if any(term in lower for term in CAP_STRATEGIES):
        return (
            "direct_5prime_capable",
            1,
            "cap-enriched 5-prime assay can test a precise initiation site",
        )

    is_long = platform in LONG_READ_PLATFORMS
    is_brain = any(term in lower for term in BRAIN_TERMS)
    is_developmental = any(term in lower for term in DEVELOPMENT_TERMS)
    is_transcriptomic = strategy in TRANSCRIPTOMIC_STRATEGIES
    is_small_rna = "mirna" in lower or "small rna" in lower

    if is_long and is_transcriptomic and not is_small_rna:
        rank = 2 if is_brain else 3
        context = "brain" if is_brain else "non-brain or mixed-tissue"
        return (
            "full_length_transcript_capable",
            rank,
            f"{context} cDNA long reads can test first-exon and splice support, "
            "but read starts are not exact capped TSS evidence",
        )
    if is_transcriptomic and not is_small_rna:
        rank = 4 if is_brain else (5 if is_developmental else 6)
        return (
            "splice_and_expression_capable",
            rank,
            "short-read RNA evidence can test first-exon junctions and "
            "expression, but not a precise TSS",
        )
    return (
        "not_useful_for_first_exon",
        9,
        "library type does not provide useful DNMT3A mRNA first-exon evidence",
    )


def parse_package(species_label: str, package: ET.Element) -> list[dict[str, str | int]]:
    experiment = package.find("EXPERIMENT")
    sample = package.find("SAMPLE")
    study = package.find("STUDY")
    if experiment is None:
        return []

    descriptor = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR")
    strategy = text(descriptor, "LIBRARY_STRATEGY")
    source = text(descriptor, "LIBRARY_SOURCE")
    selection = text(descriptor, "LIBRARY_SELECTION")
    layout_node = descriptor.find("LIBRARY_LAYOUT") if descriptor is not None else None
    layout = next(iter(layout_node), ET.Element("UNKNOWN")).tag if layout_node is not None else ""
    platform_node = experiment.find("PLATFORM")
    platform = next(iter(platform_node), ET.Element("UNKNOWN")).tag if platform_node is not None else ""
    instrument = text(platform_node, f"{platform}/INSTRUMENT_MODEL")

    sample_attributes = attributes(sample, "SAMPLE_ATTRIBUTES", "SAMPLE_ATTRIBUTE")
    sample_title = text(sample, "TITLE")
    sample_description = text(sample, "DESCRIPTION")
    study_title = text(study, "./DESCRIPTOR/STUDY_TITLE")
    experiment_title = text(experiment, "TITLE")
    tissue = sample_attributes.get(
        "tissue_type",
        sample_attributes.get("tissue", sample_attributes.get("organism part", "")),
    )
    stage = sample_attributes.get(
        "dev_stage",
        sample_attributes.get("developmental stage", sample_attributes.get("age", "")),
    )
    sex = sample_attributes.get("sex", "")
    searchable = " | ".join(
        (sample_title, sample_description, tissue, stage, study_title,
         experiment_title, strategy, selection)
    )
    evidence_class, priority_rank, evidence_limit = classify_evidence(
        strategy, selection, platform, searchable
    )

    rows: list[dict[str, str | int]] = []
    run_set = package.find("RUN_SET")
    for run in package.findall("./RUN_SET/RUN"):
        run_bases = int(run.get("total_bases") or run.get("bases") or 0)
        if not run_bases and run_set is not None and len(run_set.findall("RUN")) == 1:
            run_bases = int(run_set.get("bases") or 0)
        run_spots = int(run.get("total_spots") or run.get("spots") or 0)
        if not run_spots and run_set is not None and len(run_set.findall("RUN")) == 1:
            run_spots = int(run_set.get("spots") or 0)
        sra_file = run.find("./SRAFiles/SRAFile")
        rows.append({
            "species": species_label,
            "scientific_name": text(sample, "./SAMPLE_NAME/SCIENTIFIC_NAME"),
            "experiment_accession": experiment.get("accession", ""),
            "run_accession": run.get("accession", ""),
            "study_accession": (
                experiment.find("STUDY_REF").get("accession", "")
                if experiment.find("STUDY_REF") is not None else ""
            ),
            "bioproject": external_id(experiment.find("STUDY_REF"), "BioProject"),
            "biosample": external_id(sample, "BioSample"),
            "sample_title": sample_title,
            "tissue": tissue,
            "developmental_stage": stage,
            "sex": sex,
            "study_title": study_title,
            "library_strategy": strategy,
            "library_source": source,
            "library_selection": selection,
            "library_layout": layout,
            "platform": platform,
            "instrument": instrument,
            "spots": run_spots,
            "bases": run_bases,
            "approx_gigabases": f"{run_bases / 1e9:.3f}",
            "available": str(run.get("unavailable", "false") != "true").lower(),
            "sra_url": sra_file.get("url", "") if sra_file is not None else "",
            "brain_context": str(any(term in searchable.lower() for term in BRAIN_TERMS)).lower(),
            "developmental_context": str(
                any(term in searchable.lower() for term in DEVELOPMENT_TERMS)
            ).lower(),
            "evidence_class": evidence_class,
            "priority_rank": priority_rank,
            "evidence_limit": evidence_limit,
        })
    return rows


def parse_inventory(species_label: str, xml_path: Path) -> list[dict[str, str | int]]:
    root = ET.parse(xml_path).getroot()
    rows: list[dict[str, str | int]] = []
    for package in root.findall("EXPERIMENT_PACKAGE"):
        rows.extend(parse_package(species_label, package))
    return rows


def main() -> None:
    rows = [
        row
        for species, path in INPUTS.items()
        for row in parse_inventory(species, path)
    ]
    rows.sort(key=lambda row: (
        int(row["priority_rank"]),
        row["species"],
        -int(row["bases"]),
        row["run_accession"],
    ))
    RESULTS.mkdir(parents=True, exist_ok=True)
    output = RESULTS / "marsupial_sra_inventory.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    actionable = [
        row for row in rows
        if row["available"] == "true" and int(row["priority_rank"]) <= 4
    ]
    priority_output = RESULTS / "marsupial_sra_priority_targets.tsv"
    with priority_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(actionable)

    summary = {
        "input_xml_files": {key: str(value.relative_to(ROOT)) for key, value in INPUTS.items()},
        "runs": len(rows),
        "runs_by_species": dict(Counter(str(row["species"]) for row in rows)),
        "runs_by_evidence_class": dict(
            Counter(str(row["evidence_class"]) for row in rows)
        ),
        "available_priority_runs_rank_1_to_4": len(actionable),
        "exact_5prime_capable_runs": sum(
            row["evidence_class"] == "direct_5prime_capable" for row in rows
        ),
        "brain_long_read_runs": sum(
            row["brain_context"] == "true"
            and row["evidence_class"] == "full_length_transcript_capable"
            for row in rows
        ),
        "interpretation": (
            "Long-read cDNA can validate a candidate first exon and its splice "
            "connections, but only a cap-enriched assay can establish the exact TSS."
        ),
    }
    (RESULTS / "marsupial_sra_inventory_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
