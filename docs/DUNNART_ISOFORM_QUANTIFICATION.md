# Competitive dunnart DNMT3A transcript quantification

## Question

The exact-junction analysis showed more support for the internal
`XM_074300066.1`/`XM_074300067.1` first-junction architecture at P12 than at
P20 in dunnart neocortex. Competitive transcript quantification was added as
an independent, model-based check of that directional result.

## Reference and method

Ten RefSeq DNMT3A mRNA models were extracted from the genomic record
`Sminthopsis_crassicaudata_NC_133619.1_DNMT3A.gb`. The checksum-validated
dunnart Trinity transcriptome was used as a competitive background rather
than quantifying against the ten targets alone.

An exact 31-mer audit removed 22 obvious DNMT3A-like Trinity transcripts
having at least 100 distinct target 31-mers. The remaining 2,093,960
transcripts were retained as named Salmon decoys. This prevents reads from
being forced onto DNMT3A merely because no alternative transcript is present.
The threshold, target audit, and excluded-contig audit are explicit outputs.

Salmon 2.3.4 selective alignment used paired, inward-facing, stranded
libraries (`ISR`), concordant mappings only, and sequence- and GC-bias
correction. The first validation run used automatic library detection and
resolved to `ISR`; the remaining runs were locked to that type.

The transcript families used for interpretation were:

- full-length first-junction family: `XM_074300059.1`,
  `XM_074300062.1`, and `XM_074300064.1`;
- internal `066/067` family: `XM_074300066.1` and `XM_074300067.1`;
- alternative internal model: `XM_074300068.1`;
- four other annotated DNMT3A models.

Models `066` and `067` are combined because the diagnostic junction is shared
and ordinary short reads cannot reliably assign every fragment between
closely related full-transcript models.

## Result

| stage | specimen | internal estimated fragments/M input | internal fraction of assigned DNMT3A | internal/full model ratio |
|---|---:|---:|---:|---:|
| P12 | 436 | 9.8631 | 0.3915 | 3.6215 |
| P12 | 451 | 10.6921 | 0.4207 | 3.3761 |
| P12 | 458 | 15.2228 | 0.4475 | 3.8213 |
| P20 | 666 | 7.0286 | 0.3968 | 4.0355 |
| P20 | 746 | 5.0864 | 0.2944 | 2.2542 |
| P20 | 758 | 5.7065 | 0.2858 | 1.7767 |

Pooled P12 had 11.7962 estimated internal-family fragments per million input
fragments, compared with 5.8912 at P20, a 2.00-fold difference. The internal
family represented 42.15% of assigned DNMT3A fragments at P12 and 31.62% at
P20, a 1.33-fold difference. Each P12 animal exceeded each P20 animal for the
absolute rate (exact one-sided animal-label permutation P = 0.05). The
internal fractions overlap between stages (P = 0.10).

The internal/full-length family ratio did not separate the stages: the P12
replicate range was 3.38–3.82 and the P20 range was 1.78–4.04
(permutation P = 0.20). Thus Salmon
supports a higher absolute and proportional internal-model signal at P12,
but it provides weaker evidence for a specific internal-versus-full-length
switch. The exact internal/full-length first-junction ratio is also
non-significant (permutation P = 0.40), despite the stronger exact
internal/common-junction result.

## Interpretation and limits

This is a useful sensitivity analysis, not a population-level test. It agrees
directionally with the exact-junction result while showing that the magnitude
depends on the measurement definition:

- exact internal-junction pairs per million: 3.40-fold P12/P20;
- exact internal/common-junction ratio: 2.42-fold;
- Salmon internal-family fragments per million: 2.00-fold;
- Salmon internal fraction of assigned DNMT3A: 1.33-fold.

The attenuation is expected because transcript-level allocation shares
ambiguous fragments among closely related models. Only three animals are
available per stage. In addition, the very low reported target mapping
percentage is intentional in this single-locus analysis: tens of millions of
fragments map to the transcriptome decoys, whereas only DNMT3A target
fragments enter the reported target abundance.

No single specimen reverses the pooled Salmon direction. Across six
single-animal deletions, the internal-rate fold ranges from 1.74 to 2.18 and
the internal-fraction fold from 1.22 to 1.46. Across nine balanced
one-per-stage deletions, the ranges are 1.66–2.34 and 1.17–1.51,
respectively. The pooled internal/full-length fold also remains above one in
every deletion scenario, but individual replicate ratios still overlap
between stages; these are different claims.

Neither Salmon nor exact splice-junction counts identify a capped
transcription-start nucleotide. The data validate developmental use of an
internal splice architecture, not exact DNMT3A2 promoter initiation.

## Reproduction

Prepare the competitive reference and build its index:

```bash
python scripts/prepare_dunnart_salmon_reference.py
salmon index \
  -t data/processed/dunnart_salmon_reference/dnmt3a_refseq_targets.fasta \
     data/processed/dunnart_salmon_reference/trinity_non_dnmt3a_decoys.fasta \
  -d data/processed/dunnart_salmon_reference/decoys.txt \
  -i data/processed/dunnart_salmon_reference/salmon_index \
  -k 31 -p 8 --ramLimit 8
```

For each run, quantify with:

```bash
salmon quant \
  -i data/processed/dunnart_salmon_reference/salmon_index \
  -l ISR \
  -1 data/raw/marsupial_transcriptomes/RUN_1.fastq.gz \
  -2 data/raw/marsupial_transcriptomes/RUN_2.fastq.gz \
  -p 6 --discardOrphansQuasi --seqBias --gcBias \
  -o results/06_isoform_evolution/salmon/RUN
```

Then summarize:

```bash
python scripts/summarize_dunnart_salmon_isoforms.py
```

The principal outputs are:

- `dunnart_salmon_model_estimates.tsv`;
- `dunnart_salmon_family_estimates.tsv`;
- `dunnart_salmon_replicate_metrics.tsv`;
- `dunnart_salmon_isoform_summary.json`.
