# Mammalian DNMT3A production-tree QC

## Production analysis

- Input: 253 mammalian DNMT3A representatives, 901 codons (2,703 nt).
- Model: `GTR+F+G4` on nucleotide characters of the codon alignment.
- Constraint: monophyly of the 16 multi-species NCBI orders; seven
  singleton-sampled orders remain root-level unconstrained tips.
- Support: 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates.
- Best log likelihood: -48,716.277.
- Bootstrap split-frequency correlation: 0.999.
- Runtime: approximately 17 minutes.

The final tree contains all 253 expected tips, has nonnegative branch lengths,
and satisfies all 23 sampled order memberships in an unrooted split audit.

## Support and short branches

Among 233 internal branches not explicitly labeled as constraint nodes:

- 175 have UFBoot support at least 95%;
- 40 have UFBoot support from 70% through 94%;
- 18 have UFBoot support below 70%.

IQ-TREE reports 24 near-zero internal branches. Thirteen terminal branches are
at IQ-TREE's minimum reported length of `1e-6`; this includes the two exact
sequence-identity groups but also several tips with no uniquely assignable
DNMT3A changes. Branch-specific claims are not identifiable on these branches.

SH-aLRT values of 100% on constraint-enforced order branches are structural:
IQ-TREE warns that both alternative NNIs violate the constraint. They must not
be presented as independent sequence support.

## Composition and long-branch cautions

The IQ-TREE chi-square composition screen flags 23 of 253 sequences. These
lineages were retained because automatic deletion would disproportionately
remove marsupials and several other phylogenetically informative taxa.
Composition-aware or reduced-taxon sensitivity analyses are required before
interpreting affected long branches.

The longest terminal branches include *Erinaceus europaeus*,
*Muscardinus avellanarius*, *Puma concolor*, *Suncus etruscus*, and
*Carlito syrichta*. Their sequence models and local alignment should be checked
before any branch-selection result is prioritized.

## Diagnostic run not used

An unconstrained `GY+F3X4+G4 --fast` diagnostic was stopped during ML-distance
calculation. Its codon-model parameters reached boundaries and its provisional
topology broke five mammalian orders. It is retained only as an audit trail in
`results/03_tree/mammals_diagnostic/`.

For HyPhy, internal IQ-TREE support and constraint labels are removed from a
copy of the final tree without changing tip labels, topology, or branch
lengths. The labeled IQ-TREE tree remains the support-reporting source.

Exact sequence duplicates are then reduced to one lexicographically selected
representative per identity group for selection-model fitting. This removes
three redundant copies and yields 250 informative tips. The full 253-tip tree
and all species mappings are retained; deduplication changes no observed codon
pattern and is recorded in `results/04_selection/mammals/deduplication_audit.tsv`.

For ancestral-state interpretation only, the 250-tip tree is rerooted on the
split separating the two sampled monotremes from therians. The rooting audit
confirms that all unrooted splits and pairwise tip distances are preserved.
Site-level likelihood results are unchanged, while reconstructed change
directions now have a biologically meaningful mammalian root.
