#!/usr/bin/env python3
"""Extract total DNMT3A-like expression from the source-study Trinity counts."""

from __future__ import annotations

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from prepare_dunnart_salmon_reference import (
    automaton_shared_kmer_count,
    kmer_automaton,
    kmers,
    target_records,
)
from summarize_dunnart_neocortex_stages import exact_permutation_p_greater


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/marsupial_transcriptomes"
RESULTS = ROOT / "results/06_isoform_evolution"
ASSEMBLY = RAW / "GSE161274_FTD_Trinity_assembly.fasta.gz"
TRANSCRIPT_MAP = RAW / "GSE161274_FTD_Trinity_gene_transcript_map.txt.gz"
COUNTS = RAW / "GSE161274_FTD_gene_raw_counts.matrix.txt.gz"
JUNCTIONS = RESULTS / "dunnart_neocortex_stage_junctions.tsv"
THRESHOLD = 100
SENSITIVITY_THRESHOLDS = (100, 250, 500, 1000)
SAMPLES = {
    "FTD_P12_436.genes.results": ("SRR13036043", "P12", "436"),
    "FTD_P12_451.genes.results": ("SRR13036044", "P12", "451"),
    "FTD_P12_458.genes.results": ("SRR13036045", "P12", "458"),
    "FTD_P20_666.genes.results": ("SRR13036046", "P20", "666"),
    "FTD_P20_746.genes.results": ("SRR13036047", "P20", "746"),
    "FTD_P20_758.genes.results": ("SRR13036048", "P20", "758"),
}


def clean_header(value: str) -> str:
    return value.strip().strip('"')


def cpm(count: float, library_total: float) -> float:
    return count * 1_000_000 / library_total


def main() -> None:
    target_kmers = set().union(*(kmers(str(record.seq)) for record in target_records()))
    automaton = kmer_automaton(target_kmers)
    candidate_transcripts = {}
    with gzip.open(ASSEMBLY, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            hits = automaton_shared_kmer_count(str(record.seq), automaton)
            if hits >= THRESHOLD:
                candidate_transcripts[record.id] = {
                    "transcript_id": record.id,
                    "length": len(record.seq),
                    "distinct_dnmt3a_target_31mer_hits": hits,
                }

    transcript_to_gene = {}
    gene_to_transcripts = defaultdict(list)
    with gzip.open(TRANSCRIPT_MAP, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for gene_id, transcript_id in reader:
            transcript_to_gene[transcript_id] = gene_id
            gene_to_transcripts[gene_id].append(transcript_id)

    audit_rows = []
    candidate_genes = set()
    for transcript_id, row in candidate_transcripts.items():
        gene_id = transcript_to_gene.get(transcript_id, "")
        if gene_id:
            candidate_genes.add(gene_id)
        audit_rows.append({
            **row,
            "source_study_gene_id": gene_id,
            "mapped_to_count_matrix_gene": bool(gene_id),
        })
    audit_rows.sort(
        key=lambda row: (-row["distinct_dnmt3a_target_31mer_hits"], row["transcript_id"])
    )
    if not candidate_genes:
        raise SystemExit("No exact-sequence DNMT3A candidates map to source-study genes")

    library_totals = defaultdict(float)
    candidate_counts = defaultdict(float)
    per_gene_counts = {}
    with gzip.open(COUNTS, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = [clean_header(value) for value in next(reader)]
        sample_columns = header[1:]
        if set(sample_columns) != set(SAMPLES):
            raise SystemExit(f"Unexpected count-matrix columns: {sample_columns}")
        for values in reader:
            gene_id = clean_header(values[0])
            counts = [float(value) for value in values[1:]]
            for sample, value in zip(sample_columns, counts):
                library_totals[sample] += value
                if gene_id in candidate_genes:
                    candidate_counts[sample] += value
            if gene_id in candidate_genes:
                per_gene_counts[gene_id] = dict(zip(sample_columns, counts))

    expression_rows = []
    for sample in sample_columns:
        run, stage, specimen = SAMPLES[sample]
        expression_rows.append({
            "run_accession": run,
            "stage": stage,
            "specimen": specimen,
            "source_matrix_column": sample,
            "source_assigned_gene_count_total": round(library_totals[sample], 6),
            "dnmt3a_candidate_gene_count": round(candidate_counts[sample], 6),
            "dnmt3a_candidate_gene_cpm": round(
                cpm(candidate_counts[sample], library_totals[sample]), 6
            ),
        })

    with JUNCTIONS.open() as handle:
        junction_rows = {
            row["run_accession"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    for row in expression_rows:
        junction = junction_rows[row["run_accession"]]
        internal_rate = float(
            junction["internal_067_first_junction_pairs_per_million"]
        )
        row["internal_junction_pairs_per_million"] = internal_rate
        row["internal_junction_rate_per_dnmt3a_gene_cpm"] = round(
            internal_rate / row["dnmt3a_candidate_gene_cpm"], 9
        )

    stage_summary = {}
    for stage in ("P12", "P20"):
        rows = [row for row in expression_rows if row["stage"] == stage]
        total_counts = sum(row["dnmt3a_candidate_gene_count"] for row in rows)
        total_library = sum(row["source_assigned_gene_count_total"] for row in rows)
        cpms = [row["dnmt3a_candidate_gene_cpm"] for row in rows]
        normalized_internal = [
            row["internal_junction_rate_per_dnmt3a_gene_cpm"] for row in rows
        ]
        pooled_internal_rate = (
            sum(
                float(junction_rows[row["run_accession"]][
                    "internal_067_first_junction_pairs"
                ])
                for row in rows
            )
            * 1_000_000
            / sum(
                float(junction_rows[row["run_accession"]]["paired_read_count"])
                for row in rows
            )
        )
        pooled_gene_cpm = cpm(total_counts, total_library)
        stage_summary[stage] = {
            "replicate_count": len(rows),
            "pooled_cpm": round(pooled_gene_cpm, 6),
            "replicate_cpm_range": [round(min(cpms), 6), round(max(cpms), 6)],
            "replicate_cpms": cpms,
            "pooled_internal_junction_rate_per_gene_cpm": round(
                pooled_internal_rate / pooled_gene_cpm, 9
            ),
            "replicate_internal_junction_rate_per_gene_cpm_range": [
                round(min(normalized_internal), 9),
                round(max(normalized_internal), 9),
            ],
            "replicate_internal_junction_rates_per_gene_cpm": normalized_internal,
        }
    fold = stage_summary["P12"]["pooled_cpm"] / stage_summary["P20"]["pooled_cpm"]
    permutation_p = exact_permutation_p_greater(
        stage_summary["P12"]["replicate_cpms"],
        stage_summary["P20"]["replicate_cpms"],
    )
    normalized_fold = (
        stage_summary["P12"]["pooled_internal_junction_rate_per_gene_cpm"]
        / stage_summary["P20"]["pooled_internal_junction_rate_per_gene_cpm"]
    )
    normalized_permutation_p = exact_permutation_p_greater(
        stage_summary["P12"]["replicate_internal_junction_rates_per_gene_cpm"],
        stage_summary["P20"]["replicate_internal_junction_rates_per_gene_cpm"],
    )
    threshold_sensitivity = {}
    for threshold in SENSITIVITY_THRESHOLDS:
        genes = {
            transcript_to_gene[transcript_id]
            for transcript_id, row in candidate_transcripts.items()
            if row["distinct_dnmt3a_target_31mer_hits"] >= threshold
            and transcript_id in transcript_to_gene
        }
        sample_cpms = {
            sample: cpm(
                sum(per_gene_counts[gene][sample] for gene in genes),
                library_totals[sample],
            )
            for sample in sample_columns
        }
        p12 = [sample_cpms[sample] for sample in sample_columns if SAMPLES[sample][1] == "P12"]
        p20 = [sample_cpms[sample] for sample in sample_columns if SAMPLES[sample][1] == "P20"]
        pooled = {}
        for stage in ("P12", "P20"):
            stage_samples = [
                sample for sample in sample_columns if SAMPLES[sample][1] == stage
            ]
            pooled[stage] = cpm(
                sum(
                    per_gene_counts[gene][sample]
                    for gene in genes for sample in stage_samples
                ),
                sum(library_totals[sample] for sample in stage_samples),
            )
        threshold_sensitivity[str(threshold)] = {
            "candidate_gene_count": len(genes),
            "P12_over_P20_pooled_cpm_fold": round(
                pooled["P12"] / pooled["P20"], 6
            ),
            "all_P12_cpms_exceed_all_P20_cpms": min(p12) > max(p20),
            "animal_level_exact_permutation_p_greater": (
                exact_permutation_p_greater(p12, p20)
            ),
        }
    summary = {
        "source": "GSE161274 source-study RSEM/Trinity gene raw-count matrix",
        "candidate_rule": (
            f"Original-assembly transcript has at least {THRESHOLD} distinct exact "
            "31-mers shared with one or more dunnart RefSeq DNMT3A models."
        ),
        "candidate_transcript_count": len(candidate_transcripts),
        "candidate_count_matrix_gene_count": len(candidate_genes),
        "candidate_gene_ids": sorted(candidate_genes),
        "stages": stage_summary,
        "P12_over_P20_pooled_cpm_fold": round(fold, 6),
        "all_P12_cpms_exceed_all_P20_cpms": (
            min(stage_summary["P12"]["replicate_cpms"])
            > max(stage_summary["P20"]["replicate_cpms"])
        ),
        "animal_level_exact_permutation_p_greater": permutation_p,
        "internal_junction_rate_normalized_to_source_gene_cpm": {
            "P12_over_P20_pooled_fold": round(normalized_fold, 6),
            "all_P12_values_exceed_all_P20_values": (
                min(stage_summary["P12"][
                    "replicate_internal_junction_rates_per_gene_cpm"
                ])
                > max(stage_summary["P20"][
                    "replicate_internal_junction_rates_per_gene_cpm"
                ])
            ),
            "animal_level_exact_permutation_p_greater": normalized_permutation_p,
        },
        "exact_31mer_threshold_sensitivity": threshold_sensitivity,
        "interpretation": (
            "This source-study gene-level count is an independent measure of "
            "total DNMT3A-like expression, not an isoform or promoter measure."
        ),
    }

    with (RESULTS / "dunnart_source_dnmt3a_candidate_audit.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=audit_rows[0].keys())
        writer.writeheader()
        writer.writerows(audit_rows)
    with (RESULTS / "dunnart_source_dnmt3a_expression.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=expression_rows[0].keys()
        )
        writer.writeheader()
        writer.writerows(expression_rows)
    (
        RESULTS / "dunnart_source_dnmt3a_expression_summary.json"
    ).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
