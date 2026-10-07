# DNMT3L archive audit and DNMT3A adaptations

Reference inspected on 2026-07-22:

an external DNMT3L methodological archive that is not distributed here

## Patterns retained

- Primary inference uses a standardized, codon-aware MACSE-refined alignment.
- For the 475-sequence DNMT3A panel, a MAFFT protein guide is back-translated
  from exact-matched CDS and passed to MACSE `refineAlignment`. This avoids
  MACSE's quadratic de novo distance initialization while preserving codon-aware
  refinement; the pre-refinement projection is retained for sensitivity.
- The full-panel refinement uses one conservative leaf-cut iteration focused on
  local regions (`optim=1`, `local_realign_init=0.1`). Exhaustive branch-cut
  refinement was rejected as computationally disproportionate for 475 highly
  conserved full-length sequences; alignment sensitivity remains mandatory.
- Every chosen protein is paired to a CDS whose translation matches exactly.
- Human residue coordinates are mapped through the alignment rather than
  assumed to equal alignment columns.
- IQ-TREE and HyPhy inputs, commands, logs, and parsed tables are retained.
- The production topology uses the archive-tested `GTR+F+G4` nucleotide model
  on the codon alignment. Full codon ModelFinder is optional: in the DNMT3L
  archive it required roughly five days for 321 taxa, and the DNMT3A diagnostic
  codon fit reached parameter boundaries. HyPhy still performs codon-model
  selection inference on the chosen topology.
- FUBAR, MEME, and aBSREL are interpreted jointly, not as interchangeable tests.
- Trait and methylome analyses are explicitly associative.
- Raw public data and provenance are separated from curated intermediates.

## Changes made for DNMT3A

- DNMT3A1-like and DNMT3A2-like products are not pooled. The coding pipeline
  initially selects a full-length representative so the N-terminal region can
  be analyzed; isoform/promoter evolution will be a separate module.
- Mammalian DNMT3A is separated from teleost `dnmt3aa`/`dnmt3ab` paralogy.
- A minimum 85% length gate is applied relative to human DNMT3A1 to exclude
  DNMT3A2-like and truncated models from the full-length coding analysis; all
  exclusions are auditable.
- Domain-level summaries and a human-coordinate map are first-class outputs.
- Brain-development phenotypes are kept downstream of sequence analysis and
  require phylogenetic models, coverage thresholds, and sensitivity analyses.
- No DNMT3L residue-specific structural or life-history hypothesis was copied.

## Important archive limitation avoided

The archive's simple raw rerun used a protein alignment followed by direct
codon threading. The mature archived analysis correctly promoted standardized
MACSE output to the primary coordinate system. This project likewise treats
the MACSE-refined result—not the initial protein-guided projection—as primary.
