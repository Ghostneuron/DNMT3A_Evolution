# DNMT3A project status and inference boundaries

## Overall status

The computational discovery phase is complete. Ortholog curation, codon-aware
alignment, tree inference, four-alignment sensitivity analysis, corrected
full-site MEME, annotation-aware masking, candidate ranking, lineage
reconstruction, developmental-expression integration, human cell-type mCH
context, functional-variant prioritization, and regulatory-tail subregion
profiling are implemented and tested. In-silico functional prioritization now
also includes a 22-metric disordered-tail property screen and masked-marginal
scoring with two ESM-2 model sizes. Promoter analysis now covers both the
DNMT3A2-like internal promoter and the upstream DNMT3A1 full-length-transcript
TSS cluster. A 253-species isoform-distribution audit and a 470-way
great-ape-consensus human-branch promoter screen are also complete.
An isoform-functional module now audits 231 public GEO samples from three
matched perturbation studies. Reanalysis of 13,159,458 shared CpGs in P21
cortical neurons and 296,070 MM285 probes across E15.5 and P21 supports a
developmental division of labor between DNMT3A2 and DNMT3A1.
An independent all-context reanalysis of eight P21 NeuN-positive EM-seq
libraries now directly measures mCA: `Dnmt3a1` knockout reduces corrected
autosomal mCA by approximately 99.4%, `Dnmt3a2` knockout retains bulk mCA, and
the DNMT3A1 N-terminal deletion retains approximately 40% of WT. Gene-body
effects were integrated with all 16 matched whole-cortex RNA count tables.
Independent 2025 neuron-subclass gene sets show that core and recurrent
MeCP2-repressed genes are strongly enriched among the largest absolute
DNMT3A1-dependent mCA losses, primarily because those genes carry high WT mCA.
Four checksum-validated P18 cortex ChIP-seq tracks were also integrated with
the P21 methylomes. WT DNMT3A occupancy predicts mCA preservation after
DNMT3A2 loss, while H3K27me3-rich gene bodies show a small but chromosome-
robust excess loss after DNMT3A1 N-terminal deletion, including after joint
adjustment for WT DNMT3A occupancy.
The previously unmodeled P21 NeuN-positive H2AK119ub track has now been added.
It confirms that N-terminal-deletion mCA loss is concentrated in Polycomb-rich
genes, but does not distinguish H2AK119ub from correlated H3K27me3. A direct
gene-wise occupancy-change-to-mCA-change relationship is not detected.
A benchmarked protein-interface screen has now mapped naturally observed UDR
substitutions onto two independent cryo-EM models (PDB 8U5H and 8QZM). All
published loss-of-binding controls were directionally recovered. L188H was
concordantly predicted to strengthen the protein interface and P186Q to weaken
it, while M185V was template-discordant. Site-saturation mutagenesis supports
the L188H and P186Q directions but also shows template-sensitive magnitudes.
Native H2A-K119–ubiquitin-G76 WT, P186Q, and L188H simulation systems have now
been reconstructed from PDB 8U5H using published Amber LYX/GLX parameters.
All three explicit-solvent systems passed clash minimization and restrained
2 ps 50–300 K topology sanity checks. Independent buffer solvation produced
different water counts and box dimensions, so their total energies are not a
mutation comparison. All three also passed a 2 ps PME/NPT protocol smoke test
through complete restraint release while retaining the native isopeptide
link. Each system has now completed a single 100 ps unrestrained PME/NPT
relaxation with ten sparse frames. All three maintain the linkage and stable
temperature/density. Candidate contact differences are descriptive only and
are not converged mutation-effect estimates. Two additional 50 ps runs are
complete for all three systems: one reseeded from the pilot endpoints and one
started from the shared preproduction state with fresh velocities. Four of six
candidate/partner distance directions agree across all three runs, whereas no
candidate has an all-three-run-concordant UDR RMSD shift. H188 forms
side-chain hydrogen bonds in every run, but its partner changes. Measured CPU
throughput is approximately 1.1--2.35 ns/day, depending on the system.
A subsequent WT--L188H branch with variant-specific solvent boxes adds 500 ps
per state. H188 is farther from H2A and closer to ubiquitin in both halves of
that branch. An independently prepared common-solvent-box 500 ps pair does not
replicate those directions: H188 is instead consistently closer to H3, while
H2A and ubiquitin differences are small. UDR RMSD and residue-188 RMSF also
reverse direction between preparations. The cross-preparation structural
hypothesis is therefore sampling of alternative local interface states, not a
reproducible partner-distance, mobility, or binding effect.

The manuscript focus is now evolutionary constraint within the DNMT3A1
regulatory tail. The experimentally supported 164–219
nucleosome/H2AK119ub-engagement region is strongly conserved relative to the
variable upstream tail. G34, T12, and S97 remain secondary hypotheses for
functional follow-up, not corrected alignment-independent discoveries. None
of the six tier-1 substitutions is exceptional after within-variant correction
in the property screen, and none is unusually incompatible in either ESM-2
model relative to recurrent natural alternatives.

## What can be claimed now

- DNMT3A contains sharply different evolutionary regimes across its domains
  and within the DNMT3A1-specific N-terminal tail.
- Only 6 of 56 retained sites in residues 164–219 are polymorphic, compared
  with 132 of 154 sites in residues 1–163 (`p = 3.26e-24`).
- G34, T12, and S97 are the most defensible secondary coding candidates after
  alignment, annotation, and occupancy sensitivity analysis.
- Adult human primary-motor-cortex neurons have approximately fourfold higher
  median mCH than non-neuronal nuclei in both available donors.
- Developing dunnart neocortex shows a replicated directional increase in the
  internal DNMT3A transcript signal at P12 relative to P20.
- Human and mouse have directly CAGE-supported, syntenically corresponding
  DNMT3A1 TSS clusters. Mapped annotated starts approach within 18 bp, and one
  mapped human dominant CTSS exactly matches a mouse dominant CTSS.
- Downstream-start-compatible products occur in 13 of 18 sampled placental
  orders and in the marsupial orders Dasyuromorphia and Diprotodontia. This is
  compatible with a therian-ancestral internal-transcript architecture.
- Six recurrent lineage-supported substitutions define a prespecified
  hypothesis panel, but current sequence-only tests do not prioritize one as a
  functional perturbation.
- In GSE295720 brain, `Dnmt3a2` loss has the larger methylation effect at E15.5
  (-0.02265 mean beta), whereas `Dnmt3a1` loss has the larger effect at P21
  (-0.07031). A matched P21 cortical-neuron EM-seq analysis independently
  estimates -0.07836 for `Dnmt3a1` knockout and +0.00355 for `Dnmt3a2`
  knockout across shared CpGs.
- In the all-context P21 neuronal reanalysis, DNMT3A1 is required for nearly
  all bulk mCA above lambda background, while the DNMT3A1 N-terminal deletion
  has a reproducible partial-loss phenotype.
- Greater DNMT3A1-dependent gene-body mCA loss has a very weak inverse
  association with matched whole-cortex expression change; its adjusted
  partial R-squared is below 0.002 and does not support strong genome-wide
  transcriptional coupling.
- Core and recurrent MeCP2-repressed genes are 3.63- and 10.11-fold enriched,
  respectively, in the strongest decile of DNMT3A1-knockout mCA loss.
- Among genes with WT corrected mCA at least 0.005, WT DNMT3A-FLAG enrichment
  explains 11.5% of adjusted residual mCA-effect variation after `Dnmt3a2`
  loss. In the DNMT3A1 N-terminal-deletion model, the joint DNMT3A/H3K27me3
  chromatin features explain 1.10% of adjusted residual variation, with
  H3K27me3 associated with greater loss.
- In precisely age- and cell-matched P21 neuronal chromatin, H2AK119ub alone
  predicts greater N-terminal-deletion mCA loss (beta -0.0394), but its joint
  coefficient is effectively zero after H3K27me3 adjustment. H3K27me3 remains
  negative (beta -0.0558).
- Two independent structure templates prioritize L188H and P186Q as plausible
  interface-tuning substitutions; this is a protein-interface hypothesis, not
  a demonstrated methylation or developmental phenotype.
- In single 100 ps relaxation pilots, both variants retain the UDR--histone
  interface. L188H is locally WT-like at H2A/H3, while P186Q shows shorter
  nearest H2A/ubiquitin distances despite its static destabilization score.
  This method disagreement supports replication rather than a directional
  binding claim.
- Pilot residue-pair mapping nominates H2A N89/E91 and ubiquitin K48 around
  Q186, and H2A E91/N94, H3 R134/A135, and ubiquitin A46/G47 around H188.
  These are monitoring targets for replicas, not validated altered contacts.
- Histone-aligned pilot UDR RMSD is higher for P186Q (1.69 Å) than WT
  (1.00 Å) or L188H (1.19 Å), but neither candidate has a concordant RMSD
  direction across all three runs. Mobility differences are therefore not a
  reproduced mutation effect.
- Across all three runs, Q186 remains closer to H2A and ubiquitin, while H188
  remains WT-like at H2A and closer to ubiquitin. H3 distance directions are
  variable, and the P186/Q186--H3 separation is noncontacting.
- H188 makes side-chain hydrogen bonds in every run, but its partner switches
  from H3 A135 in the first two runs to H2A E91 (with a sparse H3 R134 contact)
  in the third. This supports altered local polar-contact capacity, not a
  unique persistent bond. Coupled with the static screen, L188H remains the
  cleaner structural hypothesis without establishing altered affinity.
- In the variant-specific-box 500 ps extension, H188 is farther from H2A than WT L188 by
  0.642 Å and closer to ubiquitin by 0.250 Å. Both directions occur in each
  250 ps half; H3 direction varies. The longer result supersedes the short-run
  description of H188 as H2A-WT-like and supports local interface
  redistribution rather than uniformly stronger binding.
- In the independent common-solvent-box 500 ps comparison, H188 is closer to
  H3 by 1.331 Å overall, with the same direction in each half (-1.115 and
  -1.547 Å). H2A (-0.102 Å) and ubiquitin (-0.031 Å) shifts are small. H3 R134
  and A135 are within 4.5 Å of H188 in all 20 frames, but no persistent H188
  side-chain hydrogen bond is detected.
- No material partner-distance or mobility direction is reproduced across
  both 500 ps preparations. This preparation sensitivity supports an
  ensemble-shifting hypothesis, not a fixed structural mechanism.

## What cannot be claimed

- No candidate is demonstrated to alter DNMT3A activity or regulation.
- No candidate is associated with brain mCH or developmental expression across
  species.
- No result demonstrates brain adaptation.
- Sequence-property and protein-language-model scores do not measure DNMT3A
  activity, phosphorylation, chromatin binding, or methylation.
- The dunnart data do not establish an internal-over-full-length isoform switch
  or an exact capped marsupial transcription start.
- Promoter sequence similarity does not yet establish orthologous promoter
  function across Theria.
- Canonical DNMT3A1 promoter-cluster detection in brain samples does not
  establish brain specificity or enrichment.
- Neither DNMT3A promoter is enriched for fixed human-branch substitutions
  relative to aligned local flanks in the strict four-great-ape screen.
- Annotation non-detection cannot be interpreted as lineage-specific isoform
  loss.
- The isoform-functional results now measure mCH, but they do not connect
  either isoform to brain size, IQ, or evolutionary change.
- The weak RNA--mCA association is cross-assay and not cell- or animal-matched;
  it does not establish that mCA loss caused the expression changes.
- MeCP2-target enrichment in absolute mCA loss does not demonstrate
  preferential fractional targeting: WT mCA is already much higher in these
  genes, and thresholded retention sensitivity tests do not support a larger
  fractional loss caused by the DNMT3A1 N-terminal deletion.
- Cortex chromatin associations do not demonstrate direct recognition of
  H3K27me3, identify a causal recruitment mechanism, or provide animal-level
  replication. The ChIP tracks are single source-study profiles from P18,
  whereas methylation was measured at P21.
- The matched P21 H2AK119ub analysis does not rescue a direct-recognition
  claim: H2AK119ub is highly correlated with H3K27me3, contributes no
  independent N-tail-deletion mCA coefficient in the joint model, and measured
  DNMT3A occupancy change does not predict mCA change gene by gene.
- EvoEF2 scores are not physical binding free energies and omit DNA, solvent,
  and the nonstandard H2A–ubiquitin linkage. They cannot establish altered
  chromatin occupancy, catalytic activity, mCA, or brain adaptation.
- The WT, P186Q, and L188H MD sanity checks establish topology stability only.
  They are not equilibrated ensembles, binding free energies, or evidence that
  either substitution changes nucleosome recognition. Their independently
  solvated boxes also prevent comparison of absolute total energies.
- The 100 ps contact summaries are not equilibrated binding ensembles. Ten
  frames from one trajectory per variant cannot establish altered affinity,
  contact occupancy, chromatin recruitment, methylation, or adaptation.
- The reseeded 50 ps runs share pilot-1 endpoint coordinates. They test
  sensitivity to velocities but are not independently prepared structural
  replicas and do not supply inferential replication.
- The third run restarts from the common protocol-defined preproduction state,
  but it remains only 50 ps and the variants retain separately solvated boxes.
  The three-run summaries are sensitivity evidence, not converged ensembles or
  statistically independent estimates of binding energetics.
- The 500 ps WT--L188H extensions increase timescale coverage but remain one
  branch per state, inherit replicate-3 endpoint preparation, and use
  variant-specific solvent boxes. Their RMSD and contact differences are
  descriptive, not replicated mutation effects or binding free energies.
- The common-solvent 500 ps comparison controls box composition and starting
  solvent coordinates but remains one short trajectory per state. Agreement
  between its two halves is not independent replication, and disagreement
  with the earlier 500 ps preparation precludes a directional binding or
  mobility claim.

## Highest-information next actions

1. Obtain biological replicates or raw-read-level replicate information for
   WT and N-tail-deleted neuronal DNMT3A occupancy. The available single-track
   comparison cannot support animal-level recruitment inference.
2. Acquire comparable neuronal mCH or developmental-expression measurements
   from replicated non-reference lineages, especially *Acomys*,
   stenodermatine bats, and monotremes.
3. Test the prespecified young-TE association after obtaining the original
   Zoonomia table S4 workbook; do not substitute figure-derived values.
4. Generate cap-enriched marsupial 5-prime data and perform an orthology-aware
   therian promoter analysis.
5. Test whether the developmental DNMT3A2-to-DNMT3A1 methylation handoff is
   conserved across mammals using matched cell types and developmental stages.
6. If additional computation is justified, extend multiple independently
   prepared WT and L188H common-box replicas beyond 500 ps or use a validated
   alchemical free-energy design. Current sampling totals 2.6 ns, but the two
   500 ps preparations disagree on the dominant local geometry and do not
   provide a converged mutation-effect estimate.

Additional post hoc selection tests on the same alignments should not change
candidate status without new homology, functional, or phenotype evidence.

Revised manuscript framing and the new regional analysis:

- `docs/MANUSCRIPT_FOCUS.md`
- `docs/REGULATORY_TAIL_EVOLUTION.md`
- `docs/IN_SILICO_FUNCTIONAL_TESTS.md`
- `docs/TE_ASSOCIATION_PLAN.md`
- `docs/DNMT3A1_PROMOTER_EVOLUTION.md`
- `docs/DENSE_PROMOTER_ISOFORM_EVOLUTION.md`
- `docs/DNMT3A_ISOFORM_FUNCTIONAL_COMPARISON.md`
- `results/04_selection/mammals/regulatory_tail/summary.json`

Machine-readable status:

- `results/08_synthesis/project_readiness.tsv`
- `results/08_synthesis/project_readiness.json`
