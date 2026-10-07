# Evidence for therian ancestry of the DNMT3A internal promoter

This document concerns the internal promoter associated with DNMT3A2-like
transcripts. The upstream full-length DNMT3A1 promoter has now been analyzed
separately in `docs/DNMT3A1_PROMOTER_EVOLUTION.md`.

## Current conclusion

The most parsimonious working hypothesis is now that a DNMT3A internal
transcript architecture was present before the eutherian–marsupial split.
This hypothesis is substantially better supported than it was at the initial
protein-only screen, but it is not yet a demonstrated ancestral functional
promoter.

## Evidence layers

### 1. Conserved coding-entry neighborhood

Human, mouse, tammar wallaby, and fat-tailed dunnart internal transcripts all
enter the full-length DNMT3A coding structure within the narrow amino-acid
neighborhood 212–222. This was determined from genomic CDS-exon overlap, not
from protein names.

The transcript structures differ in detail:

- human has two target-specific coding segments before the shared anchor;
- mouse and tammar enter directly into a shared coding segment after an
  upstream untranslated exon;
- the two dunnart models use one or two transcript-specific exons before the
  shared anchor.

These differences are compatible with lineage-specific remodeling around a
conserved internal regulatory architecture.

### 2. Direct eutherian initiation evidence

Human and mouse internal starts are supported by FANTOM5 CAGE. Whole-genome
chain mapping places the complete human internal-promoter window across the
mouse internal TSS, and one human dominant CTSS maps exactly to the mouse
annotated TSS and dominant CTSS.

Thus, the human–mouse internal regulatory interval is supported by transcript
annotation, direct capped-RNA initiation, and whole-genome synteny.

### 3. Marsupial RNA-supported transcript models

The marsupial `XM_` records are computational Gnomon models, but they are not
protein-only predictions:

- tammar `XM_072633891.1`: one supporting mRNA and 43 long SRA reads;
- dunnart `XM_074300067.1`: 148 long SRA reads, 100% RNA-seq feature coverage,
  and all annotated introns supported in 23 samples;
- dunnart `XM_074300068.1`: 138 long SRA reads and 98% RNA-seq feature
  coverage.

All three also carry poly(A) coordinate evidence. These qualifiers support the
transcript bodies and splice structures, but poly(A) evidence concerns the
3′ end and does not establish the exact capped 5′ TSS.

Targeted reanalysis of public long-read resources now adds evidence
independent of these annotation qualifiers. A dunnart multi-tissue Trinity
transcript contains five exact 35-mers spanning the predicted
`XM_074300067.1` first-exon junction. Four raw tammar ONT cDNA reads contain
exact 35-mers spanning the `XM_072633891.1` first-exon junction; one aligns
through model transcript base 1. These results validate internal first-exon
splice architecture in two marsupial orders, but the non-cap-enriched
libraries still do not establish exact transcription initiation.

Adult-cerebellum ONT data add a useful boundary condition: DNMT3A is directly
expressed in this brain region, but 20 diagnostic reads support a full-length
first junction and none supports either predicted internal first junction.
This tissue-specific negative does not weaken the independent marsupial
internal-transcript evidence from tammar gonad and the dunnart multi-tissue
assembly, but it prevents claiming constitutive or adult-cerebellum usage.

Developing-neocortex short-read data now supply direct brain-context evidence.
The internal first junction shared by `XM_074300066.1` and
`XM_074300067.1` is supported by nine unique read pairs in P20 replicate 1
and five in P20 replicate 2. The five exact junction seeds occur only in those
two internal models among ten annotated DNMT3A transcripts and never in the
full-length models. This establishes replicated use of an internal DNMT3A
transcript architecture in developing marsupial neocortex. It does not yet
establish the exact capped start. At that point, matched P12 libraries were
needed to test developmental-stage regulation.

The matched P12 analysis is now complete at three animals per stage.
Internal-junction rates are 0.6080–0.8626 pairs per million at P12, compared
with 0.1167–0.2777 at P20. The pooled P12/P20 contrast is 3.40-fold by
library depth and 2.42-fold after normalization to a common downstream DNMT3A
junction. Every P12 animal exceeds every P20 animal for both metrics (exact
one-sided animal-label permutation P = 0.05). This adds a replicated
developmental-neocortex association to the evidence for a functional
marsupial internal transcript architecture. It still does not convert
ordinary RNA-seq junctions into direct evidence of a capped promoter start.

Competitive quantification against ten DNMT3A models and a decoy-filtered
dunnart transcriptome independently supports the direction of this result:
the internal-`066/067` estimate per input fragment is 2.00-fold higher at
P12, and all three P12 animals exceed all three P20 animals for this absolute
rate (P = 0.05). The
internal/full-length-model ratio nevertheless overlaps between stages.
Accordingly, the model-based result is supportive, whereas exact diagnostic
junctions remain the more specific evidence for a developmental
splice-architecture difference. See
`docs/DUNNART_ISOFORM_QUANTIFICATION.md`.

### 4. Marsupial upstream-sequence conservation

The 2-kb upstream windows from tammar and dunnart share 254–275 unique 15-mers
and a complex 72-bp exact block. The two dunnart starts are 375 bp apart, so
their windows substantially overlap; they are alternative models from one
locus rather than independent species observations.

By contrast, exact human–marsupial upstream similarity is weak: only 1–9
shared unique 15-mers, with the shortest comparison dominated by a GC-rich
repeat. This does not refute deep orthology, because unconstrained promoter
sequence can erode over the therian divergence, but it means that
human–marsupial promoter-sequence orthology has not been demonstrated.

### 5. Monotreme inference remains open

Platypus contains the homologous full-length coding neighborhood but no
distinct internal-start transcript in its current RefSeq annotation. This is
an annotation absence, not a biological loss. It currently prevents placing
the origin confidently before the therian ancestor.

## Evidence-weighted interpretation

The data now support two related statements at different confidence levels:

1. **High confidence:** homologous DNMT3A internal promoters operate in human
   and mouse.
2. **Moderate-to-high confidence:** an internal transcript architecture near
   the same coding anchor existed in the therian ancestor.

The second statement remains a phylogenetic inference. Distinguishing ancestry
from convergent recruitment of the same coding neighborhood requires at least
one of:

- cap-enriched 5′-end evidence in a marsupial;
- a validated cross-therian whole-genome alignment through the internal
  regulatory interval;
- functional promoter assays or chromatin initiation marks;
- informative monotreme transcript data.

## Reproducible outputs

- `internal_start_exon_architecture.tsv`
- `promoter_exact_sequence_similarity.tsv`
- `marsupial_transcript_model_evidence.tsv`
- `marsupial_sra_inventory.tsv`
- `dunnart_transcriptome_first_exon_evidence.tsv`
- `tammar_ont_first_exon_evidence.tsv`
- `dunnart_P20_neocortex_junction_replicates.tsv`

All are under `results/06_isoform_evolution`. Rebuild with:

```bash
python scripts/internal_start_architecture.py
python scripts/promoter_sequence_similarity.py
python scripts/marsupial_model_evidence.py
python scripts/inventory_marsupial_sra.py
python scripts/dunnart_transcriptome_first_exon.py
python scripts/tammar_ont_first_exon.py
python scripts/summarize_dunnart_neocortex_junctions.py
python -m unittest discover -s tests
```

Full read-level methods and limitations are in
`docs/MARSUPIAL_5PRIME_EVIDENCE.md`.
