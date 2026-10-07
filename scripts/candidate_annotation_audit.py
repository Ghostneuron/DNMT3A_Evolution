#!/usr/bin/env python3
"""Audit exon support and local homology for disputed DNMT3A N-terminal sites."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from difflib import SequenceMatcher
from pathlib import Path

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/04_selection/mammals/candidate_annotation_validation"
PROTEINS = ROOT / "results/01_curated/DNMT3A_representative_proteins.faa"
PRIMARY_AA = ROOT / "results/02_alignment/DNMT3A_MACSE_AA.fasta"
SITES = {"P5": 5, "S6": 6, "E18": 18, "G34": 34}
RELATIVES = {
    "Puma_concolor": [
        "Acinonyx_jubatus",
        "Herpailurus_yagouaroundi",
        "Felis_catus",
        "Prionailurus_bengalensis",
    ],
    "Vombatus_ursinus": [
        "Phascolarctos_cinereus",
        "Notamacropus_eugenii",
        "Trichosurus_vulpecula",
    ],
}


def read_one_fasta(path: Path) -> str:
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"Expected one FASTA record in {path}, found {len(records)}")
    return str(records[0].seq).upper()


def exact_cds_offset(rna: str, cds: str) -> int:
    starts = [match.start() for match in re.finditer(f"(?={re.escape(cds)})", rna)]
    if len(starts) != 1:
        raise ValueError(f"Expected one exact CDS occurrence in RNA, found {len(starts)}")
    return starts[0]


def parse_exonerate(path: Path) -> dict[str, object]:
    text = path.read_text()
    gene_match = re.search(
        r"\texonerate:est2genome\tgene\t.*?identity ([0-9.]+) ; similarity ([0-9.]+)",
        text,
    )
    similarity = next(
        (
            line
            for line in text.splitlines()
            if "\texonerate:est2genome\tsimilarity\t" in line
        ),
        None,
    )
    if gene_match is None or similarity is None:
        raise ValueError(f"Could not parse exonerate GFF in {path}")

    # Exonerate's Align triplets are target_start, query_start, aligned_length.
    query_exons = [
        (int(query_start), int(query_start) + int(length) - 1)
        for _, query_start, length in re.findall(
            r"Align\s+(\d+)\s+(\d+)\s+(\d+)", similarity
        )
    ]
    donors = re.findall(r"\tsplice5\t.*?splice_site \"([A-Z]+)\"", text)
    acceptors = re.findall(r"\tsplice3\t.*?splice_site \"([A-Z]+)\"", text)
    return {
        "identity_pct": float(gene_match.group(1)),
        "similarity_pct": float(gene_match.group(2)),
        "query_exons": query_exons,
        "donors": Counter(donors),
        "acceptors": Counter(acceptors),
    }


def locate_interval(
    start: int, end: int, exons: list[tuple[int, int]]
) -> tuple[list[int], int | None]:
    exon_numbers = [
        index
        for index, (exon_start, exon_end) in enumerate(exons, start=1)
        if start <= exon_end and end >= exon_start
    ]
    if len(exon_numbers) != 1:
        return exon_numbers, None
    exon_start, exon_end = exons[exon_numbers[0] - 1]
    return exon_numbers, min(start - exon_start, exon_end - end)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def human_site_columns(alignment: dict[str, str]) -> dict[str, int]:
    human = alignment["Homo_sapiens"]
    wanted = {position: label for label, position in SITES.items()}
    columns: dict[str, int] = {}
    residue_number = 0
    for column, residue in enumerate(human):
        if residue == "-":
            continue
        residue_number += 1
        if residue_number in wanted:
            columns[wanted[residue_number]] = column
    if set(columns) != set(SITES):
        raise ValueError("Could not map every human candidate site in primary alignment")
    return columns


def site_classification(
    species: str, site: str, aligned_residue: str
) -> tuple[str, str]:
    if aligned_residue == "-":
        return (
            "alignment_gap_no_residue",
            "the production alignment has no species residue in this human-numbered column",
        )
    if species == "Puma_concolor" and site in {"P5", "S6", "E18"}:
        return (
            "exclude_from_homologous_site_inference",
            "genome-encoded predicted first exon conflicts with the conserved felid N terminus",
        )
    if species == "Puma_concolor":
        return (
            "retain_only_with_alignment_qualification",
            "second-exon residue lies inside the close-felid exact-match anchor, but upstream numbering is shifted",
        )
    return (
        "unresolved_model_homology",
        "genome-encoded predicted model has an unusual N-terminal extension and lacks independent transcript support",
    )


def main() -> None:
    manifest = {
        row["safe_id"]: row
        for row in csv.DictReader(
            (OUT / "record_manifest.tsv").open(), delimiter="\t"
        )
    }
    representative = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(PROTEINS, "fasta")
    }
    primary_alignment = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(PRIMARY_AA, "fasta")
    }
    site_columns = human_site_columns(primary_alignment)

    site_rows: list[dict[str, object]] = []
    annotation_rows: list[dict[str, object]] = []
    for species in RELATIVES:
        species_dir = OUT / species
        rna = read_one_fasta(species_dir / f"{species}.rna.fna")
        cds = read_one_fasta(species_dir / f"{species}.cds.fna")
        protein = read_one_fasta(species_dir / f"{species}.protein.faa")
        cds_offset = exact_cds_offset(rna, cds)
        exonerate = parse_exonerate(species_dir / "exonerate_est2genome.txt")
        exons = exonerate["query_exons"]
        assert isinstance(exons, list)

        annotation_rows.append({
            "species": species,
            "annotation_release_date": manifest[species]["annotation_release_date"],
            "assembly_accession": manifest[species]["assembly_accession"],
            "scaffold_type": manifest[species]["scaffold_type"],
            "transcript_type": manifest[species]["transcript_type"],
            "independent_transcript_evidence": "no",
            "rna_genome_identity_pct": exonerate["identity_pct"],
            "rna_genome_similarity_pct": exonerate["similarity_pct"],
            "exon_count": len(exons),
            "donor_sites": ";".join(
                f"{key}:{value}" for key, value in sorted(exonerate["donors"].items())
            ),
            "acceptor_sites": ";".join(
                f"{key}:{value}" for key, value in sorted(exonerate["acceptors"].items())
            ),
            "cds_start_in_rna_1based": cds_offset + 1,
            "cds_end_in_rna_1based": cds_offset + len(cds),
        })

        aligned_sequence = primary_alignment[species]
        for label in SITES:
            column = site_columns[label]
            aligned_residue = aligned_sequence[column]
            aa_position = sum(
                residue != "-" for residue in aligned_sequence[:column + 1]
            )
            decision, reason = site_classification(
                species, label, aligned_residue
            )
            if aligned_residue == "-":
                rna_start = rna_end = ""
                exon_numbers = []
                boundary_distance = None
            else:
                if protein[aa_position - 1] != aligned_residue:
                    raise ValueError(
                        f"{species} {label}: alignment/protein residue mismatch"
                    )
                rna_start = cds_offset + (aa_position - 1) * 3 + 1
                rna_end = rna_start + 2
                exon_numbers, boundary_distance = locate_interval(
                    rna_start, rna_end, exons
                )
            site_rows.append({
                "species": species,
                "human_numbered_site": label,
                "primary_alignment_column_1based": column + 1,
                "model_aa_position": "" if aligned_residue == "-" else aa_position,
                "model_residue": aligned_residue,
                "rna_nt_start_1based": rna_start,
                "rna_nt_end_1based": rna_end,
                "exon_number": ",".join(map(str, exon_numbers)),
                "crosses_exon_junction": (
                    "not_applicable"
                    if aligned_residue == "-"
                    else ("yes" if len(exon_numbers) != 1 else "no")
                ),
                "nearest_exon_boundary_nt": (
                    "" if boundary_distance is None else boundary_distance
                ),
                "genome_encoded_in_model": "no_site" if aligned_residue == "-" else "yes",
                "positional_homology_decision": decision,
                "reason": reason,
            })

    anchor_rows: list[dict[str, object]] = []
    for species, relatives in RELATIVES.items():
        target = representative[species][:200]
        for relative in relatives:
            reference = representative[relative][:200]
            match = SequenceMatcher(
                None, target, reference, autojunk=False
            ).find_longest_match()
            anchor_rows.append({
                "target_species": species,
                "relative_species": relative,
                "window_aa": 200,
                "target_match_start_1based": match.a + 1,
                "relative_match_start_1based": match.b + 1,
                "exact_match_length_aa": match.size,
                "exact_match_sequence": target[match.a:match.a + match.size],
                "target_prefix_before_anchor_aa": match.a,
                "relative_prefix_before_anchor_aa": match.b,
                "prefix_length_difference_aa": match.a - match.b,
            })

    write_tsv(OUT / "annotation_model_audit.tsv", annotation_rows)
    write_tsv(OUT / "candidate_site_exon_audit.tsv", site_rows)
    write_tsv(OUT / "n_terminal_relative_anchor_audit.tsv", anchor_rows)

    summary = {
        "scope": list(RELATIVES),
        "candidate_sites": list(SITES),
        "annotation_models": annotation_rows,
        "site_decision_counts": dict(
            Counter(row["positional_homology_decision"] for row in site_rows)
        ),
        "puma_close_relative_anchor": {
            "target_start_1based": 26,
            "relative_start_1based": 57,
            "minimum_exact_length_aa": min(
                int(row["exact_match_length_aa"])
                for row in anchor_rows
                if row["target_species"] == "Puma_concolor"
            ),
            "interpretation": (
                "Puma's unique 25-aa predicted first exon replaces or omits the "
                "conserved approximately 56-aa felid prefix; its P5/S6/E18 "
                "residues are not defensible homologues of human P5/S6/E18."
            ),
        },
        "vombatus_interpretation": (
            "The wombat sequence is genome-encoded in a predicted model, but its "
            "unusual N-terminal extension cannot be distinguished from a true "
            "lineage-specific change without independent RNA or a revised annotation."
        ),
        "inference_boundary": (
            "Perfect RNA-to-source-genome alignment validates the annotation model's "
            "sequence, not orthology of an alignment column or transcript expression."
        ),
    }
    (OUT / "candidate_annotation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
