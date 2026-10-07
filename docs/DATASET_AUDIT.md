# Production DNMT3A dataset audit

## Package

- File: `data/raw/ncbi_dataset/data/ncbi_dataset_dnmt3a_ortholog_sequence.zip`
- Retrieved: 2026-07-22
- Query context: human DNMT3A (NCBI GeneID 1788) ortholog sequences
- Embedded MD5 validation: all six packaged data/metadata files passed

## Raw contents

- Gene records: 521
- Protein records: 2,432
- CDS records: 2,432
- RNA records: 2,559
- Gene metadata records: 514

Metadata symbols were overwhelmingly DNMT3A ortholog labels: 464 `DNMT3A`,
44 `Dnmt3a`, two `dnmt3a`, and four uncharacterized `LOC` symbols linked by the
NCBI ortholog group.

## Taxonomic scope

This is not a mammal-only package. It contains sarcopterygian vertebrates,
including mammals, birds, reptiles, amphibians, and the West African lungfish.
No teleost `dnmt3aa`/`dnmt3ab` expansion was observed in the package metadata.
Results must be stratified by clade after joining a trusted taxonomy table.

## Full-length curation

The production curation requires:

- exact CDS translation to the selected protein;
- protein length between 85% and 120% of human DNMT3A1 (912 aa);
- one representative per species;
- configured human anchor NP_072046.2;
- preference for RefSeq NP records, annotated isoform 1, and then proximity to
  human DNMT3A1 length.

Outcome: 475 full-length DNMT3A1-like species representatives selected from
2,432 protein/CDS isoform records. Short DNMT3A2-like and truncated models were
excluded from this full-length coding analysis but remain documented in
`results/01_curated/representative_audit.tsv`.

The selected proteins range from 779 to 1,012 aa (median 910 aa). Six selected
records use curated RefSeq `NP_` accessions and 469 use predicted `XP_`
accessions. This is expected for a broad comparative package but makes
annotation-quality sensitivity analyses essential: candidate substitutions or
branches supported only by predicted models must be checked against the genome,
transcript evidence, alignment context, and alternative representative choices.

Of the 1,957 protein records not selected, 1,298 were alternate exact-matching
isoforms, 654 failed the full-length gate, and five lacked an exact CDS
translation match. No unmatched record was admitted to the codon alignment.
