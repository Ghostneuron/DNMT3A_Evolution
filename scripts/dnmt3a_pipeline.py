#!/usr/bin/env python3
"""Reproducible DNMT3A sequence-evolution workflow.

The script intentionally uses only the Python standard library. External
bioinformatics programs are invoked as explicit, logged commands.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import shutil
import subprocess
import sys
import tomllib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(make_record(header, chunks))
                header, chunks = line[1:], []
            elif header is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                chunks.append(line)
    if header is not None:
        records.append(make_record(header, chunks))
    return records


def make_record(header: str, chunks: list[str]) -> FastaRecord:
    return FastaRecord(header.split()[0], header, "".join(chunks).upper())


def write_fasta(rows: Iterable[tuple[str, str]], path: Path, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in rows:
            handle.write(f">{name}\n")
            for start in range(0, len(seq), width):
                handle.write(seq[start:start + width] + "\n")


def header_fields(description: str) -> dict[str, str]:
    return dict(re.findall(r"\[([^=\]]+)=([^\]]+)\]", description))


def normalize_nt(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper().replace("U", "T")


def translate(seq: str) -> str:
    nt = normalize_nt(seq)
    aa = [CODON_TABLE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3)]
    if aa and aa[-1] == "*":
        aa.pop()
    return "".join(aa)


def safe_id(species: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", species).strip("_")


def load_config(path: Path) -> tuple[dict, Path]:
    with path.open("rb") as handle:
        cfg = tomllib.load(handle)
    root = path.resolve().parent.parent
    return cfg, root


def project_path(root: Path, value: str) -> Path:
    p = Path(value)
    return p if p.is_absolute() else root / p


def paths(cfg: dict, root: Path) -> dict[str, Path]:
    return {key: project_path(root, value) for key, value in cfg["paths"].items()}


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        if not fieldnames:
            return
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def command_text(command: list[str]) -> str:
    return shlex.join(command)


def run_command(command: list[str], log_path: Path, cwd: Path, dry_run: bool = False,
                stdout_path: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as log:
        log.write("$ " + command_text(command) + "\n")
        if dry_run:
            log.write("DRY RUN\n")
            print(command_text(command))
            return
        if stdout_path:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            with stdout_path.open("w") as out:
                result = subprocess.run(command, cwd=cwd, text=True, stdout=out, stderr=log)
        else:
            result = subprocess.run(command, cwd=cwd, text=True, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Command failed ({result.returncode}); see {log_path}")


def init_project(cfg: dict, root: Path) -> None:
    p = paths(cfg, root)
    for key in ("curated_dir", "alignment_dir", "tree_dir", "selection_dir"):
        p[key].mkdir(parents=True, exist_ok=True)
    (root / "data/raw/ncbi_dataset/data").mkdir(parents=True, exist_ok=True)
    (root / "data/traits").mkdir(parents=True, exist_ok=True)
    print("Initialized project directories.")


def curate(cfg: dict, root: Path) -> None:
    p = paths(cfg, root)
    protein_path, cds_path = p["protein_fasta"], p["cds_fasta"]
    for required in (protein_path, cds_path):
        if not required.exists():
            raise FileNotFoundError(f"Missing required input: {required}")

    proteins = read_fasta(protein_path)
    cds_records = read_fasta(cds_path)
    cds_by_species: dict[str, list[tuple[FastaRecord, str]]] = defaultdict(list)
    for rec in cds_records:
        fields = header_fields(rec.description)
        species = fields.get("organism", "")
        if species:
            cds_by_species[species].append((rec, translate(rec.sequence)))

    project = cfg["project"]
    target = project["gene"].casefold()
    human_len = int(project["human_length_aa"])
    min_len = human_len * float(project["minimum_length_fraction"])
    max_len = human_len * float(project["maximum_length_fraction"])
    human_acc = project["human_protein_accession"]
    candidates: dict[str, list[dict]] = defaultdict(list)
    audit: list[dict] = []

    for protein in proteins:
        fields = header_fields(protein.description)
        species = fields.get("organism", "")
        gene = fields.get("gene", fields.get("gene_synonym", ""))
        base = {
            "species": species,
            "protein_accession": protein.identifier,
            "gene_header_value": gene,
            "isoform": fields.get("isoform", ""),
            "protein_length_aa": len(protein.sequence.rstrip("*")),
            "cds_accession": "",
            "translation_match": False,
            "status": "",
            "reason": "",
        }
        if not species:
            base.update(status="rejected", reason="missing_organism_header")
            audit.append(base)
            continue
        if gene and target not in {x.casefold() for x in re.split(r"[,;/ ]+", gene) if x}:
            base.update(status="rejected", reason="non_target_gene_header")
            audit.append(base)
            continue
        plen = len(protein.sequence.rstrip("*"))
        if not (min_len <= plen <= max_len):
            base.update(status="rejected", reason="outside_length_gate")
            audit.append(base)
            continue
        matches = [(c, aa) for c, aa in cds_by_species.get(species, []) if aa == protein.sequence.rstrip("*")]
        if not matches:
            base.update(status="rejected", reason="no_exact_cds_translation_match")
            audit.append(base)
            continue
        cds_rec = sorted(matches, key=lambda item: item[0].identifier)[0][0]
        row = dict(base)
        row.update(
            cds_accession=cds_rec.identifier,
            translation_match=True,
            status="candidate",
            reason="exact_match_within_length_gate",
            protein_sequence=protein.sequence.rstrip("*"),
            cds_sequence=normalize_nt(cds_rec.sequence),
            is_refseq_np=protein.identifier.startswith("NP_"),
            is_configured_human=protein.identifier == human_acc,
        )
        candidates[species].append(row)
        audit.append(base | {
            "cds_accession": cds_rec.identifier,
            "translation_match": True,
            "status": "candidate",
            "reason": "exact_match_within_length_gate",
        })

    selected: list[dict] = []
    selected_accessions: set[str] = set()
    for species, group in sorted(candidates.items()):
        choice = max(
            group,
            key=lambda r: (
                bool(r["is_configured_human"]),
                bool(r["is_refseq_np"]),
                str(r["isoform"]) == "1",
                -abs(int(r["protein_length_aa"]) - human_len),
                int(r["protein_length_aa"]),
                str(r["protein_accession"]),
            ),
        )
        choice["safe_id"] = safe_id(species)
        selected.append(choice)
        selected_accessions.add(choice["protein_accession"])

    for row in audit:
        if row["status"] == "candidate":
            if row["protein_accession"] in selected_accessions:
                row.update(status="selected", reason="best_exact_full_length_representative")
            else:
                row.update(status="rejected", reason="alternate_exact_isoform")

    if not selected:
        raise RuntimeError("No valid protein/CDS-matched representatives were selected")
    human_species = project["human_species"]
    if human_species not in {r["species"] for r in selected}:
        raise RuntimeError(f"Configured human anchor species was not selected: {human_species}")

    out = p["curated_dir"]
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "representative_audit.tsv", audit)
    metadata_fields = [
        "safe_id", "species", "protein_accession", "cds_accession", "isoform",
        "protein_length_aa", "translation_match",
    ]
    write_tsv(out / "representatives.tsv", selected, metadata_fields)
    write_fasta(((r["safe_id"], r["protein_sequence"]) for r in selected), out / "DNMT3A_representative_proteins.faa")
    write_fasta(((r["safe_id"], r["cds_sequence"]) for r in selected), out / "DNMT3A_representative_cds.fna")
    summary = {
        "raw_proteins": len(proteins),
        "raw_cds": len(cds_records),
        "selected_species": len(selected),
        "rejected_protein_records": sum(r["status"] == "rejected" for r in audit),
        "human_anchor": human_species,
        "human_accession_configured": human_acc,
        "taxon_scope_configured": project["taxon_scope"],
    }
    (out / "curation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def align(cfg: dict, root: Path, dry_run: bool) -> None:
    p = paths(cfg, root)
    curated, out = p["curated_dir"], p["alignment_dir"]
    cds = curated / "DNMT3A_representative_cds.fna"
    proteins = curated / "DNMT3A_representative_proteins.faa"
    if not cds.exists():
        raise FileNotFoundError(f"Run curate first; missing {cds}")
    out.mkdir(parents=True, exist_ok=True)
    log = out / "alignment_commands.log"
    if log.exists() and not dry_run:
        log.unlink()
    mafft_aa = out / "DNMT3A_MAFFT_AA_guide.fasta"
    reuse_mafft = (
        not dry_run and mafft_aa.exists() and mafft_aa.stat().st_size > 0
        and mafft_aa.stat().st_mtime >= proteins.stat().st_mtime
    )
    if reuse_mafft:
        log.write_text(f"Reusing current MAFFT guide: {mafft_aa}\n")
    else:
        run_command(
            [cfg["tools"]["mafft"], "--auto", str(proteins)],
            log, root, dry_run=dry_run, stdout_path=mafft_aa,
        )
    projected = out / "DNMT3A_MAFFT_projected_codons.fasta"
    reuse_projection = (
        not dry_run and projected.exists() and projected.stat().st_size > 0
        and projected.stat().st_mtime >= max(mafft_aa.stat().st_mtime, cds.stat().st_mtime)
    )
    if not dry_run and not reuse_projection:
        thread_cds_to_protein_alignment(mafft_aa, cds, projected)
    macse_command = [
        cfg["tools"]["macse"], "-prog", "refineAlignment", "-align", str(projected),
        "-out_NT", str(out / "DNMT3A_MACSE_NT.fasta"),
        "-out_AA", str(out / "DNMT3A_MACSE_AA.fasta"),
        "-max_refine_iter", str(cfg["alignment"]["macse_refine_iterations"]),
        "-optim", str(cfg["alignment"]["macse_optimization"]),
        "-local_realign_init", str(cfg["alignment"]["macse_local_realign_initial_fraction"]),
    ]
    run_command(macse_command, log, root, dry_run=dry_run)


def thread_cds_to_protein_alignment(protein_alignment: Path, cds_fasta: Path, output: Path) -> None:
    """Project exact-matched CDS codons through a protein alignment."""
    cds_by_id = {record.identifier: normalize_nt(record.sequence) for record in read_fasta(cds_fasta)}
    rows: list[tuple[str, str]] = []
    for protein in read_fasta(protein_alignment):
        if protein.identifier not in cds_by_id:
            raise ValueError(f"No curated CDS found for aligned protein {protein.identifier}")
        nt = cds_by_id[protein.identifier]
        codons = [nt[i:i + 3] for i in range(0, len(nt) - 2, 3)]
        if codons and CODON_TABLE.get(codons[-1]) == "*":
            codons.pop()
        projected: list[str] = []
        codon_index = 0
        for residue in protein.sequence:
            if residue == "-":
                projected.append("---")
                continue
            if codon_index >= len(codons):
                raise ValueError(f"Protein alignment consumes more codons than available for {protein.identifier}")
            codon = codons[codon_index]
            translated = CODON_TABLE.get(codon, "X")
            if residue != "X" and translated != residue:
                raise ValueError(
                    f"Codon projection mismatch for {protein.identifier} at residue {codon_index + 1}: "
                    f"{codon}->{translated}, alignment has {residue}"
                )
            projected.append(codon)
            codon_index += 1
        if codon_index != len(codons):
            raise ValueError(
                f"Protein alignment left {len(codons) - codon_index} CDS codons unused for {protein.identifier}"
            )
        rows.append((protein.identifier, "".join(projected)))
    write_fasta(rows, output)


def load_domains(root: Path) -> list[tuple[str, int, int]]:
    rows: list[tuple[str, int, int]] = []
    with (root / "config/human_domains.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append((row["domain"], int(row["start_aa"]), int(row["end_aa"])))
    return rows


def domain_for(position: int | None, domains: list[tuple[str, int, int]]) -> str:
    if position is None:
        return "unmapped"
    return next((name for name, start, end in domains if start <= position <= end), "other")


def clean(cfg: dict, root: Path) -> None:
    p = paths(cfg, root)
    source = p["alignment_dir"] / "DNMT3A_MACSE_NT.fasta"
    if not source.exists():
        raise FileNotFoundError(f"Run align first; missing {source}")
    records = read_fasta(source)
    lengths = {len(r.sequence) for r in records}
    if len(lengths) != 1:
        raise ValueError("MACSE nucleotide alignment sequences have unequal lengths")
    aln_len = lengths.pop()
    if aln_len % 3:
        raise ValueError(f"MACSE nucleotide alignment length is not divisible by 3: {aln_len}")

    clean_rows: dict[str, list[str]] = {}
    qc_rows: list[dict] = []
    for rec in records:
        codons: list[str] = []
        counts = defaultdict(int)
        for start in range(0, aln_len, 3):
            codon = rec.sequence[start:start + 3].upper().replace("U", "T")
            if codon == "---":
                codons.append(codon); counts["gap_codons"] += 1
            elif "!" in codon:
                codons.append("---"); counts["frameshift_codons_masked"] += 1
            elif "-" in codon:
                codons.append("---"); counts["partial_gap_codons_masked"] += 1
            elif any(base not in "ACGT" for base in codon):
                codons.append("---"); counts["ambiguous_codons_masked"] += 1
            elif CODON_TABLE.get(codon) == "*":
                codons.append("---"); counts["stop_codons_masked"] += 1
            else:
                codons.append(codon); counts["retained_codons"] += 1
        clean_rows[rec.identifier] = codons
        qc_rows.append({"safe_id": rec.identifier, "input_codons": len(codons), **counts})

    nseq = len(clean_rows)
    min_occ = float(cfg["project"]["minimum_clean_codon_occupancy"])
    ncols = aln_len // 3
    retained_columns = [
        col for col in range(ncols)
        if sum(clean_rows[name][col] != "---" for name in clean_rows) / nseq >= min_occ
    ]
    filtered = {name: "".join(codons[col] for col in retained_columns) for name, codons in clean_rows.items()}
    for row in qc_rows:
        retained = sum(clean_rows[row["safe_id"]][col] != "---" for col in retained_columns)
        row["retained_after_column_filter"] = retained
        row["retained_fraction_after_column_filter"] = f"{retained / len(retained_columns):.6f}"
    out = p["alignment_dir"]
    clean_path = out / "DNMT3A_MACSE_clean_HyPhy.fasta"
    write_fasta(filtered.items(), clean_path)
    write_tsv(out / "sequence_codon_qc.tsv", qc_rows)

    human_id = safe_id(cfg["project"]["human_species"])
    if human_id not in clean_rows:
        raise RuntimeError(f"Human anchor ID {human_id} missing from MACSE alignment")
    domains = load_domains(root)
    human_position = 0
    crosswalk: list[dict] = []
    retained_set = set(retained_columns)
    filtered_site = 0
    for original_col in range(ncols):
        human_codon = clean_rows[human_id][original_col]
        mapped_position: int | None = None
        if human_codon != "---":
            human_position += 1
            mapped_position = human_position
        if original_col in retained_set:
            filtered_site += 1
            occupancy = sum(clean_rows[name][original_col] != "---" for name in clean_rows) / nseq
            crosswalk.append({
                "filtered_codon_site_1based": filtered_site,
                "original_macse_codon_1based": original_col + 1,
                "human_DNMT3A1_aa": mapped_position if mapped_position is not None else "",
                "domain": domain_for(mapped_position, domains),
                "occupancy": f"{occupancy:.6f}",
                "human_codon": human_codon,
                "human_residue": CODON_TABLE.get(human_codon, "-") if human_codon != "---" else "-",
            })
    write_tsv(out / "human_coordinate_crosswalk.tsv", crosswalk)
    summary = {
        "sequences": nseq,
        "input_codon_columns": ncols,
        "retained_codon_columns": len(retained_columns),
        "minimum_occupancy": min_occ,
        "human_mapped_residues_before_column_filter": human_position,
    }
    (out / "alignment_qc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def analysis_alignment(p: dict[str, Path], dataset: str) -> Path:
    if dataset == "full":
        return p["alignment_dir"] / "DNMT3A_MACSE_clean_HyPhy.fasta"
    return p["alignment_dir"] / "subsets" / f"DNMT3A_{dataset}.fasta"


def tree_prefix(p: dict[str, Path], dataset: str) -> Path:
    label = "DNMT3A_tree" if dataset == "full" else f"DNMT3A_{dataset}_tree"
    return p["tree_dir"] / dataset / label


def tree(cfg: dict, root: Path, dry_run: bool, dataset: str) -> None:
    p = paths(cfg, root)
    alignment = analysis_alignment(p, dataset)
    if not alignment.exists():
        raise FileNotFoundError(f"Run clean and taxonomy subset construction first; missing {alignment}")
    out = tree_prefix(p, dataset).parent
    out.mkdir(parents=True, exist_ok=True)
    prefix = tree_prefix(p, dataset)
    command = [
        cfg["tools"]["iqtree"], "-s", str(alignment),
        "-st", str(cfg["tree"].get("sequence_type", "DNA")),
        "-m", str(cfg["tree"]["model"]),
        "-B", str(cfg["tree"]["bootstrap_replicates"]),
        "-alrt", str(cfg["tree"]["alrt_replicates"]),
        "-T", str(cfg["tree"]["threads"]), "--prefix", str(prefix),
    ]
    if dataset in {"mammals", "mammals_high_coverage"} and cfg["tree"].get("mammal_constraint"):
        constraint = project_path(root, str(cfg["tree"]["mammal_constraint"]))
        if not constraint.exists():
            raise FileNotFoundError(f"Build taxonomy subsets first; missing constraint {constraint}")
        command.extend(["-g", str(constraint)])
    run_command(command, out / "tree_command.log", root, dry_run=dry_run)


def selection(cfg: dict, root: Path, dry_run: bool, dataset: str) -> None:
    p = paths(cfg, root)
    alignment = analysis_alignment(p, dataset)
    hyphy_tree = Path(str(tree_prefix(p, dataset)) + "_HyPhy.nwk")
    tree_path = hyphy_tree if hyphy_tree.exists() else Path(str(tree_prefix(p, dataset)) + ".treefile")
    if dataset in {"mammals", "mammals_high_coverage"}:
        unique_alignment = p["alignment_dir"] / "subsets" / "DNMT3A_mammals_HyPhy_unique.fasta"
        unique_tree = Path(str(tree_prefix(p, "mammals")) + "_HyPhy_unique.nwk")
        rooted_unique_tree = Path(str(tree_prefix(p, "mammals")) + "_HyPhy_unique_rooted.nwk")
        if unique_alignment.exists() and unique_tree.exists():
            alignment = unique_alignment
            tree_path = rooted_unique_tree if rooted_unique_tree.exists() else unique_tree
    if not alignment.exists() or (not tree_path.exists() and not dry_run):
        raise FileNotFoundError("Run clean and tree before selection")
    out = p["selection_dir"] / dataset
    out.mkdir(parents=True, exist_ok=True)
    hyphy = cfg["tools"]["hyphy"]
    jobs: list[tuple[str, list[str]]] = []
    if cfg["selection"].get("run_fubar", True):
        jobs.append(("FUBAR", [hyphy, "fubar", "--alignment", str(alignment), "--tree", str(tree_path), "--output", str(out / "DNMT3A.FUBAR.json")]))
    if cfg["selection"].get("run_meme", True):
        candidate_sites = cfg["selection"].get("meme_candidate_sites", [])
        if candidate_sites:
            for raw_site in candidate_sites:
                site = int(raw_site)
                if site < 1:
                    raise ValueError("MEME candidate sites must be positive 1-based indices")
                # HyPhy 2.5.79's list matcher requires a comma-prefixed value.
                # Multi-digit values can also match numeric prefixes, so each
                # requested site is run and reported separately.
                jobs.append((
                    f"MEME_site_{site}",
                    [
                        hyphy, "meme", "--alignment", str(alignment), "--tree", str(tree_path),
                        "--limit-to-sites", f",{site}",
                        "--output", str(out / f"DNMT3A.site_{site}.MEME.json"),
                    ],
                ))
        else:
            jobs.append((
                "MEME",
                [
                    hyphy, "meme", "--alignment", str(alignment), "--tree", str(tree_path),
                    "--output", str(out / "DNMT3A.MEME.json"),
                ],
            ))
    if cfg["selection"].get("run_absrel", True):
        jobs.append(("aBSREL", [hyphy, "absrel", "--alignment", str(alignment), "--tree", str(tree_path), "--branches", "All", "--output", str(out / "DNMT3A.aBSREL.json")]))
    manifest = []
    for name, command in jobs:
        analysis_log = out / f"{name}.log"
        if analysis_log.exists() and not dry_run:
            analysis_log.unlink()
        if "--output" in command and not dry_run:
            output_path = Path(command[command.index("--output") + 1])
            if output_path.exists():
                output_path.unlink()
        run_command(command, analysis_log, root, dry_run=dry_run)
        manifest.append({"analysis": name, "command": command_text(command), "status": "dry_run" if dry_run else "completed"})
    write_tsv(out / "selection_manifest.tsv", manifest)


def mle_table(data: dict) -> tuple[list[str], list[list]]:
    mle = data.get("MLE", {})
    headers_raw = mle.get("headers", [])
    headers = [h[0] if isinstance(h, list) else str(h) for h in headers_raw]
    content = mle.get("content", {})
    if isinstance(content, dict):
        rows = content.get("0", next(iter(content.values()), []))
    else:
        rows = content
    return headers, rows


def find_column(headers: list[str], patterns: list[str]) -> int | None:
    lowered = [h.casefold() for h in headers]
    for pattern in patterns:
        for idx, header in enumerate(lowered):
            if pattern.casefold() in header:
                return idx
    return None


def summarize(cfg: dict, root: Path, dataset: str) -> None:
    p = paths(cfg, root)
    out = p["selection_dir"] / dataset
    crosswalk_path = p["alignment_dir"] / "human_coordinate_crosswalk.tsv"
    if not crosswalk_path.exists():
        raise FileNotFoundError(f"Missing coordinate map: {crosswalk_path}")
    with crosswalk_path.open() as handle:
        crosswalk = {int(r["filtered_codon_site_1based"]): r for r in csv.DictReader(handle, delimiter="\t")}
    site_rows: list[dict] = []
    specifications = [
        ("FUBAR", out / "DNMT3A.FUBAR.json", ["prob[alpha<beta]", "posterior probability", "prob"], [], None),
        ("MEME", out / "DNMT3A.MEME.json", [], ["p-value", "p value"], None),
    ]
    for raw_site in cfg["selection"].get("meme_candidate_sites", []):
        site = int(raw_site)
        specifications.append((
            "MEME", out / f"DNMT3A.site_{site}.MEME.json",
            [], ["p-value", "p value"], site,
        ))
    for method, json_path, posterior_patterns, p_patterns, requested_site in specifications:
        if not json_path.exists():
            continue
        data = json.loads(json_path.read_text())
        headers, rows = mle_table(data)
        alpha_idx = find_column(headers, ["alpha"])
        beta_idx = find_column(headers, ["sup>+</sup>", "beta+"]) if method == "MEME" else find_column(headers, ["beta"])
        posterior_idx = find_column(headers, posterior_patterns) if posterior_patterns else None
        p_idx = find_column(headers, p_patterns) if p_patterns else None
        for site, values in enumerate(rows, start=1):
            if requested_site is not None and site != requested_site:
                continue
            mapped = crosswalk.get(site, {})
            site_rows.append({
                "method": method,
                "filtered_codon_site_1based": site,
                "human_DNMT3A1_aa": mapped.get("human_DNMT3A1_aa", ""),
                "human_residue": mapped.get("human_residue", ""),
                "domain": mapped.get("domain", ""),
                "occupancy": mapped.get("occupancy", ""),
                "alpha": values[alpha_idx] if alpha_idx is not None and alpha_idx < len(values) else "",
                "beta": values[beta_idx] if beta_idx is not None and beta_idx < len(values) else "",
                "posterior_positive_selection": values[posterior_idx] if posterior_idx is not None and posterior_idx < len(values) else "",
                "p_value": values[p_idx] if p_idx is not None and p_idx < len(values) else "",
            })
    write_tsv(out / "site_selection_summary.tsv", site_rows, [
        "method", "filtered_codon_site_1based", "human_DNMT3A1_aa", "human_residue",
        "domain", "occupancy", "alpha", "beta", "posterior_positive_selection", "p_value",
    ])

    branch_rows: list[dict] = []
    absrel_path = out / "DNMT3A.aBSREL.json"
    if absrel_path.exists():
        data = json.loads(absrel_path.read_text())
        attrs = data.get("branch attributes", {})
        attrs = attrs.get("0", next(iter(attrs.values()), {})) if isinstance(attrs, dict) else {}
        for branch, values in attrs.items():
            branch_rows.append({
                "branch": branch,
                "corrected_p_value": values.get("Corrected P-value", ""),
                "uncorrected_p_value": values.get("Uncorrected P-value", ""),
                "LRT": values.get("LRT", ""),
                "rate_classes": json.dumps(values.get("Rate Distributions", ""), separators=(",", ":")),
            })
    write_tsv(out / "absrel_branch_summary.tsv", branch_rows, [
        "branch", "corrected_p_value", "uncorrected_p_value", "LRT", "rate_classes",
    ])
    print(f"Wrote {len(site_rows)} site rows and {len(branch_rows)} branch rows")


def doctor(cfg: dict, root: Path) -> int:
    rows = []
    ok = True
    version_args = {
        "macse": ["-version"], "mafft": ["--version"],
        "iqtree": ["--version"], "hyphy": ["--version"],
    }
    for label, executable in cfg["tools"].items():
        resolved = shutil.which(executable)
        status = "found" if resolved else "missing"
        version = ""
        if resolved:
            result = subprocess.run([resolved, *version_args.get(label, ["--version"])], text=True, capture_output=True)
            version = (result.stdout or result.stderr).strip().splitlines()[0] if (result.stdout or result.stderr).strip() else "version unavailable"
        else:
            ok = False
        rows.append({"tool": label, "configured": executable, "resolved": resolved or "", "status": status, "version": version})
    write_tsv(root / "results/tool_versions.tsv", rows)
    for row in rows:
        print(f"{row['tool']:8s} {row['status']:7s} {row['resolved']} {row['version']}")
    print(f"python   found   {sys.executable} {sys.version.split()[0]}")
    return 0 if ok else 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).resolve().parents[1] / "config/project.toml")
    parser.add_argument("--dry-run", action="store_true", help="print external commands without running them")
    parser.add_argument(
        "--dataset", default="full",
        choices=["full", "all_475", "all_high_coverage", "mammals", "mammals_high_coverage"],
        help="alignment subset for tree and selection stages",
    )
    parser.add_argument("stage", choices=["doctor", "init", "curate", "align", "clean", "tree", "selection", "summarize", "all"])
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg, root = load_config(args.config)
    if args.stage == "doctor":
        return doctor(cfg, root)
    if args.stage == "init":
        init_project(cfg, root)
        return 0
    stages = [args.stage] if args.stage != "all" else ["curate", "align", "clean", "tree", "selection", "summarize"]
    for stage in stages:
        print(f"== {stage} ==")
        if stage == "curate": curate(cfg, root)
        elif stage == "align": align(cfg, root, args.dry_run)
        elif stage == "clean": clean(cfg, root)
        elif stage == "tree": tree(cfg, root, args.dry_run, args.dataset)
        elif stage == "selection": selection(cfg, root, args.dry_run, args.dataset)
        elif stage == "summarize": summarize(cfg, root, args.dataset)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
