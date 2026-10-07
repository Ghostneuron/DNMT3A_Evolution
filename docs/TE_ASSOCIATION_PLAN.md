# Transposable-element association: status and prespecified analysis

## Rationale

Retrotransposon conflict is a published hypothesis for diversifying evolution
of the primate DNMT3A N terminus. Osmanski et al. measured total and young
transposable-element content across 248 placental mammal assemblies. Their
young-element measurements are more relevant than total genomic TE fraction
for a possible recurrent host-defense response.

## Data status

The required species-level values are reported in table S4 of:

Osmanski AB et al. (2023), *Insights into mammalian TE diversity through the
curation of 248 genome assemblies*, Science 380:eabn1430,
doi:10.1126/science.abn1430.

The exact source workbook is
`NIHMS1896914-supplement-Supplemental_Tables.xlsx`. Automated retrieval from
PMC, Europe PMC, and the publisher is unavailable because the author
manuscript is not in the text-and-data-mining open-access subset. The cited
Zenodo record (10.5281/zenodo.6498977) contains the analysis code but not the
source tables. No values have been inferred from figures or replaced with a
different TE dataset.

## Prespecified analysis

When the source workbook is available, the analysis will:

1. retain the original workbook unchanged;
2. extract species-level young total TE, LINE, SINE, LTR, DNA transposon, and
   rolling-circle proportions from table S4;
3. harmonize species names against the curated DNMT3A ortholog set and report
   every match and exclusion;
4. tabulate TE values by T12, G34, and S97 state;
5. count independent state transitions before statistical testing;
6. use phylogenetic models or transition-constrained permutation only when
   state replication is sufficient;
7. label uncorrected state-group comparisons as descriptive.

The primary TE outcomes are young total TE, young LINE, and young LTR
proportions. The other classes are secondary and require multiple-testing
correction.

## Decision rule

A candidate-site association will be treated as supportive only if it survives
phylogenetic control, outcome correction, and leave-one-clade-out sensitivity.
An association based on one species or one clade will not be interpreted as a
driving force. A null result will not exclude transient historical conflict,
because present-day young TE content may not represent the exposure at the
time of an ancestral substitution.

