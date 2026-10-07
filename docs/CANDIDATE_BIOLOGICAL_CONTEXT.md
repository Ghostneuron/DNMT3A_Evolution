# Biological context of retained coding candidates

The three retained sites—G34, T12, and S97—map to residues 1–178, annotated by
UniProt as a disordered region of human DNMT3A. All are within residues 1–213,
which are replaced in DNMT3A isoform 2, making them DNMT3A1-specific rather
than shared catalytic-protein changes.

| Site | Local human sequence | AlphaFold pLDDT | Exact UniProt PTM/variant | Nearest annotated PTM |
|---|---|---:|---|---|
| G34 | KDGEEQEEPRGKEERQEPSTT | 40.38 | none | phosphoserine S105, 71 aa |
| T12 | PAMPSSGPGDTSSSAAEREED | 43.75 | none | phosphoserine S105, 93 aa |
| S97 | LLPNGDLEKRSEPQPEEGSPA | 39.03 | none | phosphoserine S105, 8 aa |

AlphaFold confidence below 50 is very low and is consistent with the UniProt
disorder annotation. Consequently, these residues should not be interpreted
through predicted rigid contacts, binding pockets, or long-range structural
distances.

S97 is not itself an annotated phosphosite. Its proximity to experimentally
observed phosphoserine S105 provides a testable regulatory-tail hypothesis,
but does not demonstrate phosphorylation or phosphoregulation at S97.

## Experimental implications

The appropriate first functional tests are DNMT3A1-tail experiments:

- introduce representative ancestral/derived substitutions at G34, T12, and
  S97 in otherwise identical DNMT3A1 constructs;
- include DNMT3A2 or an N-terminal deletion as an isoform/domain control;
- measure protein abundance/stability, nuclear and chromatin localization,
  interaction or proximity profiles, and methylation activity;
- for S97, test whether substitutions alter the nearby S105 phosphorylation
  state or a broader tail-dependent interaction, without assuming S97 is a
  phosphosite.

Brain-development claims require neuronal or developmental assays. Generic
cell-line localization or methyltransferase activity would establish molecular
mechanism, not neural adaptation.

Current comparative phenotype coverage does not support an association test.
Only platypus among the 11 sampled non-serine S97 taxa has any available brain
phenotype in the project, and the developmental-expression panel contains only
serine at S97.

## Sources and outputs

Annotations come from the reviewed human DNMT3A
[UniProt entry Q9Y6K1](https://www.uniprot.org/uniprotkb/Q9Y6K1/entry).
Prediction confidence comes from
[AlphaFold DB Q9Y6K1](https://alphafold.ebi.ac.uk/entry/Q9Y6K1), model version
6. Retrieval provenance is recorded in
`data/provenance/dnmt3a_biological_context.json`.

- `results/04_selection/mammals/candidate_biological_context/retained_candidate_context.tsv`
- `results/04_selection/mammals/candidate_biological_context/summary.json`
