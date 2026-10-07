# Direct TSS evidence for DNMT3A internal promoters

## Question

Do the annotated internal DNMT3A transcript starts correspond to experimentally
observed capped-RNA initiation sites, and are those sites active in brain/CNS
samples?

## Data and coordinate compatibility

The audit uses the FANTOM5 phase 1+2 `fair+new` CAGE peak sets reprocessed on
GRCh38/hg38 and GRCm38/mm10. CAGE maps capped transcript 5′ ends and therefore
provides direct transcription-initiation evidence. The FANTOM5 reprocessing
study describes the peak conversion, read realignment, quality filtering, and
assembly-specific resources:

- https://doi.org/10.1038/sdata.2017.107
- https://fantom.gsc.riken.jp/5/datafiles/reprocessed/

Human primary-chromosome GRCh38 coordinates are directly compatible with the
RefSeq GRCh38.p14 target. Mouse FANTOM5 peaks were converted from GRCm38/mm10
to the RefSeq GRCm39 analysis assembly using the recorded UCSC
`mm10ToMm39.over.chain.gz` chain. Both peak interval endpoints and the dominant
CTSS were required to map uniquely and coherently.

Downloaded URLs, assemblies, retrieval dates, and SHA-256 hashes are recorded
in `data/raw/fantom5_cage/source_manifest.tsv`.

## Direct initiation result

- Human `NM_153759.3`: the closest same-strand CAGE peak begins 2 bp from the
  annotated internal TSS, with its dominant CTSS 7 bp away. Three additional
  same-strand peaks occur within 500 bp, including a high-count broader
  cluster.
- Mouse `NM_153743.5`: after GRCm38-to-GRCm39 conversion, the annotated TSS
  lies inside the closest peak and exactly matches its dominant CTSS.
- Tammar and dunnart remain annotation-only because no compatible
  species-specific direct-TSS atlas was included. They were not scored as
  negative.

These results establish that the canonical human and mouse internal DNMT3A
starts correspond to experimentally observed transcription initiation. They
do not by themselves establish that the human and mouse promoter sequences are
orthologous or that the marsupial models are expressed.

## Human–mouse regulatory synteny

The complete human 2-kb-upstream/500-bp-downstream internal-promoter window was
mapped through the UCSC GRCh38-to-GRCm38 whole-genome chain and then through
the GRCm38-to-GRCm39 chain. Both window endpoints map uniquely to chromosome
12 in the expected reversed gene orientation, and the mapped interval contains
the mouse `NM_153743.5` internal TSS.

The human annotated TSS maps 78 bp from the mouse annotated TSS. More
strikingly, one directly observed human CAGE dominant CTSS maps exactly to the
mouse annotated TSS and mouse dominant CAGE CTSS. This provides strong evidence
that the human and mouse internal transcription-initiation clusters occupy a
homologous regulatory interval. It does not imply that every base or
transcription-factor interaction is conserved.

Detailed landmark mappings are in
`results/06_isoform_evolution/human_mouse_internal_promoter_synteny.tsv`.

## Brain/CNS context

The corresponding raw FANTOM5 per-sample peak-count matrices were screened
using an explicit keyword classifier applied to sample descriptions. Raw
counts were reduced to detected/not-detected calls because library sizes were
not normalized in this module.

Both the closest internal-start peak and the local 500-bp CAGE cluster are
detected in human and mouse brain/CNS samples, including samples with
developmental keywords. This establishes brain-context activity, but not
brain-specific enrichment, developmental regulation, or isoform abundance.
Those stronger claims require normalized promoter-level expression and curated
sample ontologies.

## Reproducible commands

```bash
python scripts/direct_tss_evidence.py
python scripts/fantom5_tissue_context.py
python scripts/human_mouse_promoter_synteny.py
python -m unittest discover -s tests
```

Principal outputs:

- `results/06_isoform_evolution/direct_tss_evidence.tsv`
- `results/06_isoform_evolution/fantom5_nearby_cage_peaks.tsv`
- `results/06_isoform_evolution/fantom5_internal_promoter_sample_counts.tsv`
- `results/06_isoform_evolution/fantom5_tissue_context_summary.json`
- `results/06_isoform_evolution/human_mouse_internal_promoter_synteny.tsv`
- `results/06_isoform_evolution/human_mouse_internal_promoter_synteny_summary.json`

The next evolutionary inference gate is an orthology-aware regulatory
comparison using whole-genome alignments or validated liftOver chains,
followed by marsupial transcript-end evidence where available.
