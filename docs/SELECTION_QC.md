# Mammalian DNMT3A site-selection QC

## Analysis set

The primary selection input contains 250 unique mammalian DNMT3A coding
sequences and 901 codons. Three exact redundant sequences were removed without
removing their species from the archived 253-species tree and audit tables.
HyPhy estimates a global nonsynonymous/synonymous rate ratio of approximately
0.04, consistent with strong overall purifying selection.

FUBAR was run on the primary MACSE-refined alignment, a 227-sequence
composition-screen sensitivity set, and the independently filtered
pre-refinement MAFFT codon projection. Candidate sites use a posterior
probability threshold of 0.90. MEME was then used as a targeted episodic
selection follow-up on the union of the three FUBAR candidate sets.

## Candidate results

| Filtered site | Human DNMT3A1 | FUBAR primary | Composition pass | MAFFT projection | MEME p-value | Interpretation |
|---:|---:|---:|---:|---:|---:|---|
| 9 | T12 | 0.9807 | 0.8736 | 0.9480 | 0.000920 | Pervasive signal is composition-sensitive; targeted episodic signal remains strong |
| 30 | G34 | 0.9856 | 0.9654 | 0.9854 | 0.000663 | Only candidate supported by all three FUBAR analyses and targeted MEME |
| 107 | A114 | 0.8498 | 0.9386 | 0.8508 | 0.1246 | Composition-sensitivity-only FUBAR signal; MEME does not support it |

T12 and G34 remain significant after a conservative Bonferroni correction for
the three targeted MEME tests (`p` < 0.0033). These MEME tests are follow-up
analyses on sites selected from the same dataset, not independent replication
and not an alignment-wide discovery scan.

All three positions are in the low-complexity N-terminal regulatory region,
outside the PWWP, ADD, and catalytic domains. The human sequence contexts are:

- T12: `PAMPSSGPGDTSSSAAEREED`
- G34: `KDGEEQEEPRGKEERQEPSTT`
- A114: `GSPAGGQKGGAPAEGEGAAET`

MACSE-versus-MAFFT amino-acid agreement is 99.2% at T12 and 100% at G34 and
A114. Thus G34 is the most robust current candidate, while T12 warrants
additional composition-aware and N-terminal alignment scrutiny.

## Interpretation boundary

These results identify unusual coding evolution in the DNMT3A1-specific
N-terminal extension. They do not yet show altered DNMT3A activity, neuronal
methylation, brain development, or organismal adaptation. Those claims require
functional evidence and phylogenetically controlled developmental
transcriptomic or methylomic phenotypes.

The exact candidate table is
`results/04_selection/mammals/fubar_candidate_sensitivity.tsv`. Separate MEME
JSON files are used because HyPhy 2.5.79's site-list matcher can include shorter
numeric prefixes for multi-digit requested sites. The requested row is
extracted explicitly; unintended prefix rows are ignored.

Lineage-level ancestral-state and branch-posterior interpretation is reported
separately in `docs/CANDIDATE_LINEAGES.md`.

Domain-level profiling shows that the candidate concentration in the
N-terminal region is not significant after conditioning on amino-acid
variability. See `docs/DOMAIN_EVOLUTION.md`.

An independent PRANK codon alignment retains FUBAR posterior probabilities of
0.9531 for T12 and 0.9648 for G34. Targeted PRANK MEME p-values are 0.0449 and
0.0393, respectively; these are nominal and do not pass Bonferroni correction
for the three original PRANK candidate tests. PRANK does not support A114.
See `docs/PRANK_SENSITIVITY.md`.

A targeted post hoc aBSREL follow-up tested the two alignment-robust internal
branches carrying reconstructed G34 changes. Neither is significant:
the sampled *Rhinolophus* stem has LRT 0 and p=1, while the sampled lorisiform
stem has LRT 3.400, uncorrected p=0.06786, and Holm p=0.13572 across the two
foregrounds. This does not erase the site-level G34 result, because aBSREL and
MEME test different hypotheses, but it means there is no branch-wide
confirmation for either proposed lineage. See
`docs/CANDIDATE_LINEAGES.md`.

## Alignment-wide MEME qualification

A subsequent unfiltered MEME run tested all 901 retained codons. P5, S6, and
E18 pass 5% Benjamini-Hochberg correction on the primary alignment; P5 and E18
also pass 901-site Bonferroni correction. T12 and G34 have full-scan BH
q-values of 0.0829 and 0.0669, respectively, and therefore are not
alignment-wide discoveries at 5% FDR.

None of the three full-scan discoveries is significant under both independent
alternative alignments. P5 and S6 replicate in MAFFT but not PRANK, while E18
fails in both MAFFT and PRANK. The primary full-scan discoveries are therefore
classified as alignment-sensitive. See `docs/FULL_MEME.md`.

A conservative site-specific sensitivity masks the *Puma concolor* and
*Vombatus ursinus* codons that disagree among primary, MAFFT, and PRANK while
leaving the other 248 taxa unchanged. P5, S6, and E18 then have MEME p-values
of 0.667, 0.117, and 0.485, respectively. None is supported after ambiguity
masking, demonstrating that the primary discoveries depend on those disputed
placements.

An annotation-aware follow-up confirms that the disputed codons exist under
the selected genome-derived models but does not validate their homology.
Puma's P5/S6/E18 column states are residues 2, 3, and 15 of a unique predicted
first exon; residue 26 begins a 144-aa exact match to residue 57 in four close
felids. Those Puma states are excluded from homologous-site inference.
Wombat's unusual predicted N terminus remains unresolved without independent
RNA evidence. This supports retaining G34 and T12 only as exploratory coding
hypotheses, treating S6/P5 as lower-priority annotation-sensitive leads, and
deprioritizing E18. See `docs/CANDIDATE_ANNOTATION_VALIDATION.md`.

The same three-alignment rule was subsequently applied to all 266 retained
human-mapped N-terminal coordinates. Only 138 of 66,500 taxon-site cells
disagree, causing 113 informative primary codons to be masked. Puma, tarsier,
and wombat are the largest taxon outliers. Masked FUBAR retains T12 (0.9494)
and G34 (0.9616); A16 emerges only after masking (0.9338). Exact-site masked
MEME gives nominal p=0.0448 for T12 and p=0.0426 for G34, but both have Holm
p=0.1277 across T12/A16/G34. A16 is MEME-negative (p=0.2732). Thus T12 and
G34 remain qualified hypotheses rather than corrected discoveries. See
`docs/N_TERMINAL_RELIABILITY_MASK.md`.

The tarsier outlier has now been annotation-audited. Its 22 discrepancies are
a single human 34–59 block created when MACSE/MAFFT place unique Carlito model
residues 1–22 where PRANK places gaps. Downstream primate anchors show a
34–37-residue prefix deficit in the old predicted model, so the block is
excluded. Because Carlito was already masked at G34, the masked G34 signal does
not depend on this artifact.

Reliability-masked lineage reconstruction removes the Puma/wombat T12
terminal hypotheses and the Carlito G34 terminal hypothesis. Nine
change-coupled T12 branches and the original two internal G34 branches remain.
T12's raw branch posterior saturates 475 branches after masking and is not a
valid lineage-ranking statistic by itself. The final evidence matrix therefore
uses change-coupled branch summaries and retains only G34 and T12 for qualified
follow-up. See `docs/FINAL_CANDIDATE_EVIDENCE.md`.
