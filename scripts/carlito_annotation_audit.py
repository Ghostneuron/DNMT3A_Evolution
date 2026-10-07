#!/usr/bin/env python3
"""Audit the alignment-discordant Carlito DNMT3A N-terminal segment."""

from __future__ import annotations

import csv
import json
import sys
from difflib import SequenceMatcher
from pathlib import Path

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from candidate_annotation_audit import (  # noqa: E402
    exact_cds_offset,
    locate_interval,
    parse_exonerate,
    read_one_fasta,
    write_tsv,
)


SPECIES = "Carlito_syrichta"
OUT = ROOT / "results/04_selection/mammals/candidate_annotation_validation"
SPECIES_DIR = OUT / SPECIES
PRIMARY_AA = ROOT / "results/02_alignment/DNMT3A_MACSE_AA.fasta"
REPRESENTATIVES = ROOT / "results/01_curated/DNMT3A_representative_proteins.faa"
MASK_MANIFEST = (
    ROOT / "results/02_alignment/reliability_mask/n_terminal_mask_manifest.tsv"
)
RELATIVES = [
    "Otolemur_garnettii",
    "Nycticebus_coucang",
    "Aotus_nancymaae",
    "Callithrix_jacchus",
    "Homo_sapiens",
]
MOTIFS = ["DLEKRSEPQPEEG", "RGRLRGGLGWESS"]


def human_columns(human: str, positions: set[int]) -> dict[int, int]:
    columns: dict[int, int] = {}
    residue_number = 0
    for column, residue in enumerate(human):
        if residue == "-":
            continue
        residue_number += 1
        if residue_number in positions:
            columns[residue_number] = column
    if set(columns) != positions:
        raise ValueError("Could not map all requested human positions")
    return columns


def main() -> None:
    manifest_rows = list(
        csv.DictReader(MASK_MANIFEST.open(), delimiter="\t")
    )
    disputed = [
        row for row in manifest_rows if row["safe_id"] == SPECIES
    ]
    positions = {int(row["human_DNMT3A1_aa"]) for row in disputed}
    if len(disputed) != 22:
        raise RuntimeError(f"Expected 22 Carlito discordances, found {len(disputed)}")

    alignment = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(PRIMARY_AA, "fasta")
    }
    columns = human_columns(alignment["Homo_sapiens"], positions)
    carlito_aligned = alignment[SPECIES]

    rna = read_one_fasta(SPECIES_DIR / f"{SPECIES}.rna.fna")
    cds = read_one_fasta(SPECIES_DIR / f"{SPECIES}.cds.fna")
    protein = read_one_fasta(SPECIES_DIR / f"{SPECIES}.protein.faa")
    cds_offset = exact_cds_offset(rna, cds)
    exonerate = parse_exonerate(SPECIES_DIR / "exonerate_est2genome.txt")
    exons = exonerate["query_exons"]
    assert isinstance(exons, list)

    site_rows: list[dict[str, object]] = []
    disputed_by_position = {
        int(row["human_DNMT3A1_aa"]): row for row in disputed
    }
    for position in sorted(positions):
        column = columns[position]
        aligned_residue = carlito_aligned[column]
        model_position = sum(
            residue != "-" for residue in carlito_aligned[:column + 1]
        )
        if aligned_residue == "-" or protein[model_position - 1] != aligned_residue:
            raise RuntimeError(f"Alignment/protein mismatch at human position {position}")
        rna_start = cds_offset + (model_position - 1) * 3 + 1
        rna_end = rna_start + 2
        exon_numbers, boundary = locate_interval(rna_start, rna_end, exons)
        source = disputed_by_position[position]
        site_rows.append({
            "human_DNMT3A1_aa": position,
            "human_residue": alignment["Homo_sapiens"][column],
            "primary_alignment_column_1based": column + 1,
            "carlito_model_position": model_position,
            "carlito_residue": aligned_residue,
            "rna_nt_start_1based": rna_start,
            "rna_nt_end_1based": rna_end,
            "exon_number": ",".join(map(str, exon_numbers)),
            "crosses_exon_junction": "yes" if len(exon_numbers) != 1 else "no",
            "nearest_exon_boundary_nt": "" if boundary is None else boundary,
            "primary_codon": source["primary_codon"],
            "mafft_codon": source["mafft_codon"],
            "prank_codon": source["prank_codon"],
            "decision": "exclude_from_homologous_site_inference",
            "reason": (
                "primary/MAFFT place a shortened predicted first-exon segment "
                "against human positions where PRANK places a gap"
            ),
        })
    write_tsv(SPECIES_DIR / "discordant_site_exon_audit.tsv", site_rows)

    representatives = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(REPRESENTATIVES, "fasta")
    }
    target = representatives[SPECIES][:220]
    anchor_rows: list[dict[str, object]] = []
    for relative in RELATIVES:
        reference = representatives[relative][:220]
        match = SequenceMatcher(
            None, target, reference, autojunk=False
        ).find_longest_match()
        row: dict[str, object] = {
            "target_species": SPECIES,
            "relative_species": relative,
            "window_aa": 220,
            "target_match_start_1based": match.a + 1,
            "relative_match_start_1based": match.b + 1,
            "exact_match_length_aa": match.size,
            "exact_match_sequence": target[match.a:match.a + match.size],
            "longest_anchor_prefix_difference_aa": match.a - match.b,
        }
        for motif in MOTIFS:
            key = motif[:6]
            row[f"{key}_target_start_1based"] = target.find(motif) + 1
            row[f"{key}_relative_start_1based"] = reference.find(motif) + 1
        anchor_rows.append(row)
    write_tsv(SPECIES_DIR / "primate_anchor_audit.tsv", anchor_rows)

    record = next(
        row
        for row in csv.DictReader(
            (OUT / "record_manifest.tsv").open(), delimiter="\t"
        )
        if row["safe_id"] == SPECIES
    )
    summary = {
        "species": SPECIES,
        "annotation": {
            "assembly_accession": record["assembly_accession"],
            "annotation_release_date": record["annotation_release_date"],
            "transcript_type": record["transcript_type"],
            "scaffold_type": record["scaffold_type"],
            "independent_transcript_evidence": False,
            "rna_genome_identity_pct": exonerate["identity_pct"],
            "exon_count": len(exons),
            "cds_start_in_rna_1based": cds_offset + 1,
        },
        "discordant_positions": len(site_rows),
        "human_position_range": [min(positions), max(positions)],
        "carlito_model_position_range": [
            min(int(row["carlito_model_position"]) for row in site_rows),
            max(int(row["carlito_model_position"]) for row in site_rows),
        ],
        "all_discordant_codons_first_exon": all(
            row["exon_number"] == "1" for row in site_rows
        ),
        "all_discordant_codons_excluded": all(
            row["decision"] == "exclude_from_homologous_site_inference"
            for row in site_rows
        ),
        "conserved_anchor_offsets": {
            motif: {
                "carlito_start_1based": target.find(motif) + 1,
                "relative_start_range_1based": [
                    min(representatives[name][:220].find(motif) + 1 for name in RELATIVES),
                    max(representatives[name][:220].find(motif) + 1 for name in RELATIVES),
                ],
            }
            for motif in MOTIFS
        },
        "interpretation": (
            "Carlito's 2017 single predicted model replaces or omits roughly "
            "34-37 residues of the conserved primate DNMT3A1 prefix. MACSE and "
            "MAFFT force the model's unique residues 1-22 against human residues "
            "34-59, whereas PRANK gaps them. The codons are genome-encoded under "
            "the model but are not reliable positional homologues."
        ),
    }
    (SPECIES_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
