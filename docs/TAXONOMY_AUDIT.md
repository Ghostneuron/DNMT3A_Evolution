# NCBI taxonomy and analysis-subset audit

## Provenance and validation

The 475 selected species names were matched one-to-one to taxon IDs in the
NCBI Datasets DNMT3A gene report. Higher lineages were retrieved from NCBI
Taxonomy with EFetch on 2026-07-22 and cached at
`data/provenance/ncbi_taxonomy_efetch_2026-07-22.xml`.

All 475 taxon IDs were returned, all scientific names were represented, and
the taxonomy, alignment, and sequence-QC identifiers are identical. The
reproducible parser and subset builder is `scripts/taxonomy_subsets.py`.

## Composition

| Major clade | Representatives |
|---|---:|
| Mammalia | 253 |
| Aves | 140 |
| Squamata | 40 |
| Amphibia | 21 |
| Testudines | 18 |
| Crocodylia | 2 |
| Other sarcopterygian (lungfish) | 1 |

## Analysis subsets

- `all_475`: complete curated vertebrate panel.
- `all_high_coverage`: 473 species after the prespecified 85% sequence-coverage
  gate.
- `mammals`: 253 species.
- `mammals_high_coverage`: 253 species; no mammal fails the coverage gate.

The mammalian alignment is the primary coding-evolution dataset for later
brain-development integration. The broader panel provides ancestral context
and a taxon-sampling sensitivity analysis. A DNMT3A gene tree is not
automatically a trusted species tree; topology discordance and weakly supported
branches must be evaluated before branch-based selection tests.

Two exact coding-sequence identity groups occur in the mammalian set:
*Delphinus delphis*–*Stenella coeruleoalba*–*Tursiops truncatus* and
*Camelus dromedarius*–*Camelus ferus*. Their individual terminal substitutions
and branch-specific selection are not identifiable from DNMT3A coding sequence
alone. The groups are retained in an explicit audit table so zero-length or
collapsed tips are not mistaken for independent molecular observations.

The 253 mammals span 23 NCBI orders. A multifurcating order-level constraint is
written to `results/03_tree/constraints/NCBI_mammal_orders_constraint.nwk`.
It enforces only the well-established monophyly of each sampled mammalian order;
relationships within and among orders remain available for likelihood-based
resolution. Unconstrained and order-constrained topologies should be compared
before choosing the tree used by HyPhy.

Orders represented by one species are emitted as unconstrained root-level tips
rather than invalid unary clades. The resulting constraint contains every
mammalian alignment tip exactly once.
