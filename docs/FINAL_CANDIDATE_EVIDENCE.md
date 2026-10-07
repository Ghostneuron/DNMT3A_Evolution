# Final DNMT3A coding-candidate evidence matrix

## Role in the revised manuscript

The primary coding result is the contrast between a variable upstream
DNMT3A1-specific tail and a conserved nucleosome/H2AK119ub-engagement region at
residues 164–219. G34, T12, and S97 all lie upstream of that supported
functional interval. They are secondary experimental hypotheses rather than
the organizing result of the manuscript. See
`docs/REGULATORY_TAIL_EVOLUTION.md`.

## Current ranking

The coding analysis has passed through full-scan multiple-testing correction,
independent MAFFT, PRANK, and COBALT alignments, site-specific consensus masking,
taxon-wide N-terminal reliability masking, annotation audits, masked FUBAR and
MEME reruns, ancestral reconstruction, and targeted branch testing.

| Rank | Site | Mammal-scope full MEME BH q (909 codons) | Masked FUBAR | Masked MEME p | Adjusted masked p | Decision |
|---:|---|---:|---:|---:|---:|---|
| 1 | G34 | 0.0871 | 0.9616 | 0.0426 | 0.1277 | qualified primary coding hypothesis |
| 2 | T12 | 0.0871 | 0.9494 | 0.0448 | 0.1277 | qualified secondary coding hypothesis |
| 3 | S97 | 0.7206 | not applicable | 0.0224 full scan | 0.7206 across 909 codons | qualified mammal-scope tertiary hypothesis |
| 4 | S6 | 0.0320 | 0.8190 | 0.1167 | 0.3502 | exploratory annotation-sensitive lead |
| 5 | A16 | 0.7206 | 0.9338 | 0.2732 | 0.2732 | sensitivity-emergent; do not promote |
| 6 | P5 | 0.000396 | 0.0068 | 0.6667 | 0.9692 | artifact-dependent; deprioritize |
| 7 | E18 | 0.00758 | 0.2690 | 0.4846 | 0.9692 | artifact-dependent; deprioritize |
| 8 | A114 | 0.7206 | 0.8464 | not run | — | unsupported; stop follow-up |

The adjusted masked p-values use the relevant three-test sensitivity family.
The 909-codon BH values are from the corrected occupancy-within-mammals scan.
G34 and T12 remain hypotheses rather than corrected discoveries. P5, S6, and
E18 pass that primary-scan FDR, but their actions remain constrained by the
independent-alignment and annotation-aware sensitivity tests below.

COBALT independently supports both retained candidates: G34 has FUBAR
posterior 0.9634 and MEME p=2.45e-6 (Holm p=1.96e-5 across eight COBALT
follow-ups), while T12 has FUBAR posterior 0.9517 and MEME p=0.00402
(Holm p=0.0241). These results strengthen—not replace—the masked analysis.

## G34

G34 is supported by FUBAR in the primary, MAFFT, PRANK, and taxon-wide masked
alignments. The tarsier artifact is masked and the signal persists. The masked
MEME result is nominal (p=0.0426) but not significant after the three-follow-up
Holm correction (p=0.1277).

The masked reconstruction retains the sampled *Rhinolophus* stem G-to-A and
lorisiform stem G-to-S changes. Both already failed targeted post hoc aBSREL:
p=1.0 for the *Rhinolophus* stem and Holm p=0.1357 for the lorisiform stem.
Four change-coupled branches cross the descriptive posterior threshold only
after masking: the Pecora stem, hippopotamus, *Rhinolophus ferrumequinum*, and
Tasmanian devil. These are post hoc leads and should not trigger an expanding
series of branch tests without a prespecified biological hypothesis.

## T12

T12 retains FUBAR support in all alignments and after the reliability mask.
Masked MEME is nominal (p=0.0448; Holm p=0.1277). The Puma and wombat terminal
changes disappear after masking, while nine change-coupled branches remain,
including the eutherian, cetacean, *Acomys*, and sampled didelphine stems and
several terminal changes.

The raw masked MEME positive-class posterior exceeds 0.90 on 475 branches
because the fitted positive-rate class has nearly unit mixture weight. Raw
branch posterior is therefore not lineage-discriminating for T12. The lineage
table uses a stricter descriptive filter requiring a reconstructed
nonsynonymous change and an informative terminal state or an internal branch.
No T12 branch has formal branch-wide confirmation.

## Stopped or exploratory candidates

- S6 retains moderate FUBAR support but loses MEME support when the disputed
  Puma and wombat placements are removed.
- A16 crosses the FUBAR threshold only after masking and lacks MEME support.
- P5 and E18 depend on non-homologous or unresolved predicted N-terminal
  placements. E18 is significant under COBALT but remains deprioritized because
  its homology/annotation problem is upstream of the selection test.
- A114 lacks FUBAR or MEME support despite perfect four-alignment agreement.

S97 is now third-ranked. Its original absence from MACSE and MAFFT was caused
by filtering occupancy across 475 vertebrates before extracting mammals.
Mammal-scope MACSE, MAFFT, PRANK, and COBALT all give 100% occupancy, identical
states across 250 taxa, and FUBAR posterior 0.911–0.921. The corrected
alignment-wide MACSE scan gives p=0.02245 but BH q=0.7206 across 909 codons;
the other exact-site MEME fits give p=0.021–0.024. These fits reuse the same
site data and are not independent replications; S97 remains a qualified
hypothesis rather than a corrected discovery.

The S97 reconstruction infers serine at the mammalian root and eight
nonsynonymous changes coupled to MEME positive-class posterior ≥0.90,
including the sampled Stenodermatinae, *Acomys*, and monotreme stems. These
lineages are descriptive candidates only; no branch-wide selection test has
confirmed them.

## Biological boundary

None of these results demonstrates altered DNMT3A enzymatic activity,
DNMT3A1-specific regulation, neuronal mCH, developmental timing, or brain
adaptation. G34, T12, and S97 are suitable for targeted functional or comparative
follow-up only if that next analysis has an independent mechanistic or
phenotypic hypothesis.

G34, T12, and S97 all lie in the UniProt-annotated disordered,
DNMT3A1-specific N-terminal tail and have AlphaFold pLDDT below 50. S97 is
eight residues from annotated phosphoserine S105, but is not itself an
annotated phosphosite. Functional work should therefore test tail-dependent
regulation rather than predicted rigid structural contacts.

## Machine-readable outputs

- `results/04_selection/mammals/final_candidate_evidence.tsv`
- `results/04_selection/mammals/final_candidate_evidence.json`
- `results/04_selection/mammals/mammal_scope/full_meme/summary.json`
- `results/04_selection/mammals/mammal_scope/full_meme/all_sites.tsv`
- `results/04_selection/mammals/reliability_mask/lineages/branch_evidence.tsv`
- `results/04_selection/mammals/reliability_mask/lineages/substitution_events.tsv`
- `results/04_selection/mammals/reliability_mask/lineages/prioritized_change_coupled_branches.tsv`
- `results/04_selection/mammals/reliability_mask/lineages/high_posterior_comparison.tsv`
- `results/04_selection/mammals/reliability_mask/lineages/summary.json`
- `results/04_selection/mammals/cobalt/cobalt_candidate_sensitivity.tsv`
- `results/04_selection/mammals/cobalt/cobalt_sensitivity_summary.json`
- `results/04_selection/mammals/s97_homology_audit/summary.json`
- `results/04_selection/mammals/s97_homology_audit/lineages/summary.json`
