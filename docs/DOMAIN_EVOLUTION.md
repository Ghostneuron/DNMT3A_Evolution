# DNMT3A domain-level evolutionary profile

## Purpose

All three sensitivity-defined FUBAR candidates map to the operational
DNMT3A1 N-terminal regulatory region. This analysis tests whether that pattern
is stronger than expected after accounting for the region's size and unusually
high background variability.

The profile uses 250 unique mammalian coding sequences and 901 retained codons.
A site is called polymorphic when at least two non-gap amino-acid states are
observed. A recurrently variable site has at least two non-major observations,
reducing the influence of singleton sequence or annotation errors.

## Domain profile

| Region | Retained sites | Polymorphic | Recurrently variable | Mean normalized amino-acid entropy | Primary FUBAR candidates |
|---|---:|---:|---:|---:|---:|
| N-terminal regulatory | 266 | 161 | 132 | 0.0528 | 2 |
| PWWP extended | 150 | 42 | 16 | 0.0106 | 0 |
| PWWP–ADD linker | 48 | 9 | 8 | 0.0072 | 0 |
| ADD | 139 | 16 | 4 | 0.0041 | 0 |
| ADD–methyltransferase linker | 19 | 2 | 0 | 0.0009 | 0 |
| DNA methyltransferase | 279 | 45 | 6 | 0.0021 | 0 |

The N-terminal region contains 60.5% polymorphic sites, compared with 18.0% in
the rest of DNMT3A1 (odds ratio 7.01; one-sided hypergeometric
`p = 3.0e-35`). For recurrent variability, the fractions are 49.6% versus
5.4% (odds ratio 17.41; `p = 5.4e-51`).

Despite this variability, 237 of 266 retained N-terminal sites have primary
FUBAR posterior probability at least 0.90 for negative selection. The region is
therefore not generally unconstrained: a minority of positions account for
much of its amino-acid diversity.

## Candidate enrichment

For the primary FUBAR scan, both candidates are N-terminal. Relative to all 901
retained sites, the one-sided enrichment probability is 0.0869. Once restricted
to the 275 polymorphic sites—the sites with direct opportunity to produce a
selection signal—the probability is 0.342. Restriction to 166 recurrently
variable sites gives 0.631.

The union of all three sensitivity scans contains T12, G34, and A114:

| Site universe | N-terminal candidates | Total candidates | Enrichment p |
|---|---:|---:|---:|
| All retained sites | 3 | 3 | 0.0255 |
| Polymorphic sites | 3 | 3 | 0.199 |
| Recurrently variable sites | 3 | 3 | 0.500 |

The raw concentration is therefore explained by the N terminus containing a
disproportionate share of variable sites. This does not invalidate T12 or G34,
which have site-level evidence, but it argues against claiming that adaptive
evolution preferentially targets the N-terminal region as a whole.

PRANK initially added two FUBAR candidates, A16 and S97. A16 remains
alignment-sensitive and lacks PRANK MEME support. S97's original absence from
the primary and MAFFT-filtered coordinate systems was a broad-vertebrate
filtering artifact; mammal-scope refiltering validates it across all four
alignments.

## Biological interpretation

The structured PWWP, ADD, and catalytic regions are much more conserved than
the N terminus, consistent with strong constraint on chromatin reading,
allosteric regulation, and catalysis. The DNMT3A1 N terminus is the principal
reservoir of mammalian coding diversity and is consequently the most plausible
region for lineage-specific regulatory changes.

The appropriate hypothesis is now narrower: determine whether T12, G34, S97, or
their surrounding low-complexity sequence affects DNMT3A1-specific regulation,
stability, localization, protein interactions, or developmental expression.
The present data do not demonstrate brain adaptation or a DNMT3A1-versus-
DNMT3A2 functional difference.

Machine-readable outputs:

- `results/04_selection/mammals/domain_evolution_summary.tsv`
- `results/04_selection/mammals/site_domain_metrics.tsv`
- `results/04_selection/mammals/candidate_enrichment_tests.tsv`
- `results/04_selection/mammals/n_terminal_variability_tests.tsv`
- `results/04_selection/mammals/top_primary_fubar_sites.tsv`
