# Comparative analysis of the canonical DNMT3A1 promoter region

## Scope

The DNMT3A1 promoter is not represented by one conserved transcription start
site (TSS). Human and mouse each have multiple curated full-length transcripts
whose starts form a local cluster. The analysis therefore separates two
questions:

1. Is the broader full-length DNMT3A1 promoter region syntenically conserved?
2. Are individual TSSs and local promoter sequences conserved?

The primary sequence comparison uses one full-length transcript variant per
species. The human and mouse synteny analysis includes every curated
full-length-transcript TSS in the local cluster. Promoter windows contain
2,000 bp upstream and 500 bp downstream of the selected TSS.

## Direct initiation evidence

FANTOM5 CAGE supports both human full-length-transcript start regions:

- Human `NM_022552.5`: the annotated TSS lies within a same-strand CAGE peak
  and exactly matches its dominant CTSS.
- Human `NM_175629.2`: the closest same-strand peak interval is 70 bp away,
  and its dominant CTSS is 76 bp away.

The mouse cluster is also CAGE-supported:

- Mouse `NM_007872.5` and `NM_001271753.2` share an annotated TSS that exactly
  matches the dominant CTSS.
- Mouse `NM_001421857.1` starts within a same-strand peak, with the dominant
  CTSS 16 bp away.

Thus, full-length DNMT3A1 transcription initiates from multiple local starts
in both species. The annotated transcript-variant-1 start is not the only
appropriate landmark for comparing promoter evolution.

## Human–mouse synteny

Whole-genome chain mapping supports correspondence between the human and mouse
TSS clusters:

- both curated human TSS positions map uniquely into the mouse locus;
- the closest mapped human annotated TSS is 18 bp from a curated mouse TSS;
- the closest mapped human CAGE dominant CTSS is 1 bp from a curated mouse
  TSS;
- one mapped human CAGE dominant CTSS exactly matches a mouse dominant CTSS.

One endpoint of the 2.5-kb human variant-1 promoter window is not uniquely
mapped. Consequently, the complete arbitrary window cannot be described as a
single conserved block. The TSS-cluster landmarks provide stronger evidence
than the window endpoints for human–mouse regulatory-region correspondence.

## Exact sequence similarity

The selected human and mouse upstream windows share 32 unique 15-mers and a
31-bp longest exact block. The tammar wallaby and fat-tailed dunnart windows
share 496 unique 15-mers and an 89-bp longest exact block.

Exact human-to-marsupial and human-to-platypus similarity is weak at this
resolution: no shared unique upstream 15-mers were found in the selected
pairwise comparisons. This does not establish independent promoter origins.
Short promoter motifs can turn over while locus position and transcriptional
output remain conserved. Conversely, exact local matches do not establish
conserved regulatory effects.

The non-model mammal TSSs are based on predicted RefSeq transcript models.
Species-specific cap-enriched data are not available for tammar, dunnart, or
platypus in this analysis.

## Brain context

The union of local canonical DNMT3A1 CAGE peaks is detected in 98 of 100
keyword-classified human brain/CNS samples and 129 of 132 mouse brain/CNS
samples. Detection is also common in nonbrain samples: 1,598 of 1,729 human
samples and 920 of 941 mouse samples.

These raw within-library detection counts support broad promoter activity,
including brain activity. They do not establish brain enrichment,
developmental regulation, or comparable expression magnitude between species.

## Evidence-weighted conclusion

Human and mouse have syntenically corresponding, directly active DNMT3A1 TSS
clusters. Individual annotated starts differ, so promoter evolution should be
analyzed at the cluster level rather than by forcing transcript-variant-1 TSSs
into one-to-one correspondence.

The selected marsupial promoter windows are strongly conserved between tammar
and dunnart. Deeper promoter ancestry remains unresolved because their exact
TSSs are predicted and cross-therian exact sequence conservation is weak.

## Reproducible outputs

- `scripts/canonical_promoter_evolution.py`
- `scripts/canonical_promoter_tissue_context.py`
- `results/06_isoform_evolution/canonical_dnmt3a1_promoter/summary.json`
- `results/06_isoform_evolution/canonical_dnmt3a1_promoter/canonical_dnmt3a1_tss_cluster.tsv`
- `results/06_isoform_evolution/canonical_dnmt3a1_promoter/canonical_direct_tss_evidence.tsv`
- `results/06_isoform_evolution/canonical_dnmt3a1_promoter/human_mouse_canonical_promoter_synteny.tsv`
- `results/06_isoform_evolution/canonical_dnmt3a1_promoter/canonical_promoter_sequence_similarity.tsv`
- `results/06_isoform_evolution/canonical_dnmt3a1_promoter/canonical_fantom5_tissue_context_summary.json`

