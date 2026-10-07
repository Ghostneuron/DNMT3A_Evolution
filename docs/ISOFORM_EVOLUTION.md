# DNMT3A isoform-evolution screen

## Scope

This module performs an operational annotation screen, not promoter
reconstruction. Each mammalian protein product in the saved NCBI ortholog
package is compared with that species' selected full-length DNMT3A1-like
representative.

A product is called `DNMT3A2_like_downstream_start` when it has:

- a unique 8–45-aa leader;
- a longest exact shared block beginning at 18–28% of the full product;
- a shared block covering at least 70% of the full product and 90% of the
  shorter product.

These rules are calibrated against human NP_715640.2 (DNMT3A2) and distinguish
it from exact downstream suffixes without a unique leader and from other
truncated or alternatively spliced products.

Leaderless exact downstream suffixes are retained separately as
`leaderless_downstream_core_product`. They are compatible with an incomplete
downstream-start annotation but can also be simple truncation models; they are
not counted as strict A2-like products.

## Interpretation

An annotated A2-like product is useful evidence for a downstream-start
transcript model, especially when it is a curated RefSeq `NP_` product.
Predicted `XP_` products are kept as a separate evidence tier.

Failure to detect an A2-like product is **not** a gene or isoform loss. Genome
annotation depth, transcript sampling, and model propagation vary strongly
among species. Biological gain/loss claims require conserved downstream
promoter sequence, transcription-start evidence, splice-junction support, and
synteny-aware genome inspection.

## Reproducible command

```bash
python scripts/classify_dnmt3a_isoforms.py
python -m unittest discover -s tests
```

Machine-readable outputs are in `results/06_isoform_evolution`.

## Initial result and annotation audit

The screen detects operational A2-like products in 74 of 253 mammalian species:
one curated human RefSeq case and predicted products in 73 species. Candidate
shared-core junctions are tightly clustered at full-product positions 208–213,
which supports consistency of the sequence rule.

Detection is not a prevalence estimate. Species with detected products have
deeper protein annotation than nondetected species, and all detected candidates
in this package are eutherian. No metatherian or prototherian candidate is
recovered under the strict leader rule. Tammar wallaby does contain a predicted
691-aa leaderless downstream-core product (XP_072489992.1), making it a
high-priority transcript-start validation target rather than evidence of
absence. Fat-tailed dunnart has one leaderless and one extended-leader
downstream-core model. All three non-eutherian products map exactly to packaged
predicted `XM_` CDS/RNA models, but none has experimental TSS support.
Consequently, the current screen cannot test the ancestral mammalian origin or
infer lineage-specific loss.

`isoform_annotation_strata.tsv` reports detection by annotation-depth bin,
mammalian group, and order. `isoform_annotation_audit_summary.json` records
junction consistency and the explicit gain/loss inference gate.

The downstream-core protein models are also linked back to packaged CDS and RNA
records in `isoform_transcript_model_support.tsv`. An exact CDS translation and
a modeled 5′ UTR establish internal consistency of an NCBI transcript model,
but an `XM_`/`XP_` pair remains predicted evidence; it does not establish an
experimentally observed transcription start site.

## Targeted genomic transcript-start audit

Five gene-centered RefSeq genomic windows were inspected with strand-aware
coordinates: human, mouse, tammar wallaby, fat-tailed dunnart, and platypus.
The audit parses every annotated DNMT3A mRNA, links compatible CDS features,
records the first transcribed exon and genomic TSS, and measures its distance
from the outermost annotated DNMT3A start. It also exports 2-kb-upstream/500-bp-
downstream sequence windows in transcriptional orientation. These are
candidate promoter windows, not experimentally validated promoters.

The principal findings are:

- Human `NM_153759.3`/`NP_715640.2` has a curated internal start 90,275 bp
  downstream of the outermost annotated DNMT3A TSS. This is the canonical
  DNMT3A2 reference architecture.
- Mouse `NM_153743.5`/`NP_714965.1` has a curated internal start 90,259 bp
  downstream. The protein-only heuristic labels it leaderless rather than
  strict A2-like, illustrating why genomic transcript architecture must take
  precedence over the operational leader rule.
- Tammar `XM_072633891.1` starts 114,523 bp downstream and has a distinct
  230-nt first exon. Dunnart `XM_074300067.1` and `XM_074300068.1` start
  113,605 and 113,980 bp downstream and have distinct 202-nt and 517-nt first
  exons. These are coherent internal transcript architectures, but all are
  predicted RefSeq models.
- The current platypus annotation contains three transcripts sharing the
  outermost first exon and no downstream-core target. This is an absence of
  annotation, not evidence for biological absence.

Thus, internal DNMT3A transcript architecture is not restricted to eutherian
annotations: it is also predicted in two marsupial lineages. This raises the
priority of testing an origin before the therian ancestor, but does not yet
establish when a functional DNMT3A2 promoter arose. The next inference gate is
orthology-aware comparison of the internal regulatory interval plus direct TSS
support (for example CAGE/RAMPAGE or transcript-end evidence).

Reproduce the audit with:

```bash
python scripts/promoter_tss_audit.py
python -m unittest discover -s tests
```

The genomic window manifest and SHA-256 hashes are under
`data/raw/promoter_tss`. Detailed rows, the JSON summary, and oriented sequence
windows are under `results/06_isoform_evolution`.
`primary_internal_start_promoters.bed` and
`primary_internal_start_promoters.fasta` provide the conservative five-model
input set for overlap with direct TSS data and orthology-aware conservation
analysis. A naive global alignment of these promoter windows is not treated as
evidence of orthology because short noncoding similarities, repeats, and large
evolutionary distances can be misleading.

## Direct capped-RNA evidence

FANTOM5 CAGE peaks provide direct support for the canonical internal starts in
both species with compatible data. Human `NM_153759.3` has a same-strand peak
2 bp from its annotated TSS (dominant CTSS 7 bp away), and mouse
`NM_153743.5` exactly matches the dominant CTSS after an explicit
GRCm38-to-GRCm39 coordinate conversion. Both local internal-promoter clusters
are detected in brain/CNS samples in the per-sample CAGE matrices.

This upgrades human and mouse from transcript annotation to direct
transcription-initiation evidence. It does not validate marsupial initiation,
establish brain specificity, or prove promoter-sequence orthology. Full
methods, evidence tiers, and limitations are in `docs/DIRECT_TSS_EVIDENCE.md`.

Whole-genome chain mapping further shows that the complete human internal
promoter window spans the mouse internal TSS region. One human dominant CAGE
CTSS maps exactly to the mouse annotated TSS and dominant CAGE CTSS. Thus,
human–mouse internal-promoter orthology is supported at the syntenic regulatory
interval level; the remaining open question is how far this architecture
extends toward marsupials and monotremes.

## Therian ancestry assessment

Exon-anchored comparison shows that all five primary human, mouse, tammar, and
dunnart internal models enter the full-length DNMT3A product around amino acids
212–222. The marsupial models also have species-specific RNA evidence:
43–148 long SRA reads per model, 98–100% RNA-seq feature coverage for the two
dunnart transcripts, and complete intron support across 23 samples for one
dunnart model.

Tammar and dunnart upstream windows share 254–275 unique 15-mers and a complex
72-bp exact block. Exact human–marsupial upstream similarity is much weaker,
so cross-therian promoter-sequence orthology remains unresolved.

The evidence therefore supports a therian-ancestral internal-transcript
architecture as the leading hypothesis, not as a final gain/loss conclusion.
See `docs/THERIAN_PROMOTER_ANCESTRY.md` for the evidence ladder and explicit
inference limits.

## Direct marsupial transcript-structure follow-up

Reanalysis of public long-read resources has now recovered exact
junction-spanning evidence. A dunnart assembled transcript contains five exact
35-mers across the predicted `XM_074300067.1` first junction. Four independent
tammar ONT reads span the `XM_072633891.1` first junction, and one reaches
annotated transcript base 1 in local alignment.

This upgrades marsupial internal-transcript architecture from annotation-only
support to direct cDNA support in two orders. It does not validate the exact
capped TSS, because neither dataset is cap-enriched, and it does not establish
brain-specificity. See `docs/MARSUPIAL_5PRIME_EVIDENCE.md`.

An adult female dunnart cerebellum ONT replicate was subsequently analyzed in
full. All 56 retained candidates classify as DNMT3A rather than DNMT3B, and 20
span the full-length-model first junction. No read spans either predicted
internal first junction. Thus, adult cerebellum directly supports DNMT3A
expression but not internal-isoform usage in this replicate.

That developmental-neocortex test is now positive and biologically replicated.
Exact diagnostic seeds recover nine internal-junction read pairs in P20 animal
666 and five in P20 animal 746, while each animal also has two full-length
first-junction pairs. The internal seeds match only the closely related
`XM_074300066.1`/`XM_074300067.1` model family among annotated DNMT3A
transcripts. This establishes internal-transcript use in developing marsupial
neocortex, but not brain specificity, exact promoter initiation, or an
age-dependent change.

The matched stage comparison now shows a replicated directional pattern.
Three P12 animals have 0.6080–0.8626 internal-junction pairs per million,
whereas three P20 animals have 0.1167–0.2777. Pooled internal-junction
support is 3.40-fold higher at P12, and the internal/common-downstream
junction ratio is 2.42-fold higher, showing that the contrast is not explained
only by greater overall DNMT3A detection. All P12 animals exceed all P20
animals for both metrics (exact one-sided animal-label permutation P = 0.05).
P12 and P20 correspond principally
to infragranular and supragranular neuron generation, respectively, in the
source study. This supports developmental regulation of internal-transcript
use, while three animals per stage, sparse junction counts, and lack of capped
5′ data still preclude a precise population effect or exact promoter
assignment.

The direct internal/full-length first-junction ratio does not separate the
stages (pooled fold 1.15; exact permutation P = 0.40), and its direction
reverses in some deletion analyses. Accordingly, “developmental regulation”
here means higher internal-transcript signal and a higher internal/common
DNMT3A-junction ratio at P12; it does not mean a demonstrated
internal-over-full-length isoform switch.

Competitive Salmon quantification provides a model-based sensitivity check.
Against all ten DNMT3A models and 2,093,960 transcriptome decoys, pooled
internal-`066/067` fragments per million are 2.00-fold higher at P12, and
the internal fraction of assigned DNMT3A is 1.33-fold higher. All three P12
animals exceed all three P20 animals for the absolute rate (P = 0.05), but
the internal-fraction ranges overlap (P = 0.10). The
internal/full-length transcript-family ratio overlaps between stages, so
these estimates support the direction of the exact-junction result without
strengthening the claim to a discrete isoform switch. See
`docs/DUNNART_ISOFORM_QUANTIFICATION.md`.
