# Marsupial DNMT3A internal-transcript evidence

## Result

Public long-read resources now directly support the predicted internal
DNMT3A transcript architecture in both sampled marsupial lineages.

### Fat-tailed dunnart

The checksum-validated 2024 multi-tissue Trinity assembly was searched with
35-nt seeds from the two predicted internal models.

- `TRINITY_DN6521_c0_g1_i9` contains the `XM_074300067.1` internal first
  exon and the shared downstream exon.
- Five distinct exact 35-mers span the predicted first-exon splice junction.
  This is direct assembled-cDNA support for that junction, independent of the
  later RefSeq/Gnomon model.
- For `XM_074300068.1`, `TRINITY_DN390442_c0_g2_i1` contains first-exon and
  shared-exon sequence and five exact seeds spanning the second-to-shared-exon
  junction. The first-to-second-exon junction is not exactly recovered.

The assembly combines short-read RNA-seq with PacBio Iso-Seq, so it supports
transcript structure but is not an independent raw-read count.

### Tammar wallaby

All 4,590,698 reads from ONT PCR-cDNA run `SRR30198375` were scanned in both
orientations for `XM_072633891.1` first-exon, shared-exon, and junction seeds.

- 46 reads contain at least one exact 35-mer seed.
- Four reads contain exact first-junction-spanning seeds and also contain
  first-exon and shared-exon seeds.
- Three reads contain multiple independent exact junction seeds; the fourth
  contains one.
- Local alignment places one junction-supporting read at model transcript
  base 1. This is consistent with a full-length molecule reaching the
  annotated 5′ boundary.

These observations directly validate the internal first-exon splice
architecture in an independent marsupial order.

## What this establishes

The evidence materially strengthens the inference that a DNMT3A internal
transcript architecture predates the divergence of sampled marsupial orders
and is compatible with origin in the therian ancestor. The result is no
longer based only on predicted RefSeq transcript qualifiers.

## What this does not establish

The tammar library is PCR-cDNA and the dunnart assembly is not cap-enriched.
Ordinary cDNA molecules can be 5′ truncated, and their read starts do not mark
the capped initiation nucleotide. Therefore:

- the precise marsupial TSS remains unvalidated;
- promoter activity itself is not directly measured;
- these gonadal/multi-tissue data do not establish brain expression or brain
  specificity of the internal isoform.

Direct marsupial CAGE, RAMPAGE, capped long-read sequencing, or targeted
5′-RACE would still be required for exact TSS validation. The public dunnart
adult-cerebellum ONT datasets are valuable for brain-context confirmation,
but they remain ordinary cDNA and are each approximately 26–31 Gb.

## Adult cerebellum follow-up

The smallest complete adult-cerebellum replicate, `ERR15827373`, was downloaded
and checksum-validated. It contains 24,189,458 complete female cerebellum ONT
cDNA reads. One additional header-only source record was audited and skipped.

An exact-seed scan retained 56 candidate reads. Protein-level comparison
against dunnart DNMT3A and DNMT3B references classified all 56 as DNMT3A:

- 20 reads contain exact 35-mers spanning the full-length
  `XM_074300059.1` first junction;
- 36 contain downstream DNMT3A sequence but do not reach a diagnostic first
  junction;
- zero contain an exact first-junction seed for `XM_074300067.1` or
  `XM_074300068.1`.

Therefore, DNMT3A is directly detected in adult dunnart cerebellum, but this
replicate supports the full-length transcript and does not detect either
predicted internal transcript. This is a tissue-, sex-, stage-, and
replicate-specific negative. It does not demonstrate that the internal
transcripts are absent from the developing brain, other brain regions, other
sexes, or dunnart generally.

The result makes developmental neocortex RNA-seq or capped developmental-brain
5′ data more informative next targets than additional adult-cerebellum
replicates.

## Developing neocortex follow-up

Three checksum-validated P20 neocortex RNA-seq replicates from independent
animals were scanned with exact 35-nt seeds spanning diagnostic DNMT3A splice
junctions. All three replicate the predicted internal first-junction
architecture:

- `SRR13036046` (specimen 666; 42,002,168 read pairs): nine unique pairs span
  the internal junction and two span a full-length first junction;
- `SRR13036047` (specimen 746; 42,831,584 read pairs): five unique pairs span
  the internal junction and two span a full-length first junction.
- `SRR13036048` (specimen 758; 72,017,744 read pairs): 20 unique pairs span
  the internal junction and two span a full-length first junction.

Thus, 34 of 156,851,496 pooled read pairs support the internal first junction
in three biological replicates. All five diagnostic seeds match only
`XM_074300066.1` and `XM_074300067.1` among the ten annotated dunnart DNMT3A
transcript models; none matches a full-length model. Conversely, the
full-length seed set matches only `XM_074300059.1`, `XM_074300062.1`, and
`XM_074300064.1`. This validates use of the internal first-exon splice
architecture in developing marsupial neocortex and excludes a
single-library artifact.

The shared junction cannot distinguish the downstream differences between
models `066` and `067`, and ordinary RNA-seq does not validate the exact
capped TSS. Sparse junction counts are not whole-transcript abundance
estimates. Most importantly, P20 replication by itself does not establish a
developmental change relative to the adult cerebellum, because age, brain
region, sex, and sequencing technology all differ. A matched P12 neocortex
comparison was therefore selected as the appropriate stage test.

## Replicated P12–P20 developmental pattern

The matched comparison is now complete for three biological replicates per
stage. P12 corresponds to the source study's principal period of
infragranular neuron generation, whereas P20 corresponds to supragranular
neuron generation.

| stage | specimen | read pairs | internal pairs/M | internal/common ratio |
|---|---:|---:|---:|---:|
| P12 | 436 | 105,262,158 | 0.6080 | 0.7356 |
| P12 | 451 | 103,550,608 | 0.7533 | 0.8478 |
| P12 | 458 | 92,747,882 | 0.8626 | 0.7143 |
| P20 | 666 | 42,002,168 | 0.2143 | 0.2903 |
| P20 | 746 | 42,831,584 | 0.1167 | 0.2778 |
| P20 | 758 | 72,017,744 | 0.2777 | 0.3390 |

All three P12 animals therefore exceed all three P20 animals for both
prespecified metrics. Pooled descriptively, P12 has 222 internal-junction
pairs in 301,560,648 fragments, versus 34 in 156,851,496 P20 fragments. This
is a 3.40-fold difference in internal pairs per million. The common downstream
DNMT3A junction is 1.40-fold higher at P12, but internal support rises
disproportionately: the internal/common ratio is 2.42-fold higher at P12.

An exhaustive animal-label permutation test evaluates the six animals rather
than treating read pairs as biological replicates. The one-sided exact
P value is 0.05 for both internal pairs per million and the internal/common
ratio, the minimum possible with a balanced 3+3 design and 20 label
assignments.

The direct internal/full-length first-junction comparison is substantially
weaker. Its pooled P12/P20 fold is 1.15, replicate ranges overlap
(P12 5.71–8.67; P20 2.50–10.00), and the exact one-sided permutation
P value is 0.40. The internal/full-length direction also reverses in some
single-animal and balanced deletion scenarios. Thus the data establish more
internal-transcript signal at P12, including a robust increase relative to a
common downstream DNMT3A junction, but do not establish a discrete
internal-over-full-length isoform switch.

This supports a replicated directional association between developmental
stage and internal DNMT3A transcript use in dunnart neocortex. With three
animals per stage, sparse junction counts, and a minimally resolved exact
P value, the effect size remains an estimate from a small study rather than a
definitive population parameter. A targeted junction assay in an independent
cohort is the appropriate next validation. These data still do not establish
a capped TSS or direct promoter activity.

## Single-animal robustness and metadata audit

The result does not depend on one specimen. After deleting each animal in
turn, the P12/P20 internal-junction-rate fold remains 2.89–4.46 and the
internal/common-junction fold remains 2.30–2.67. In all nine balanced
analyses that delete one P12 and one P20 animal, the corresponding ranges are
2.67–4.88 and 2.23–2.78. Every scenario remains directionally above one for
these two metrics. In contrast, the exact internal/full-length fold ranges
from 0.90 to 1.87 under single-animal deletion and from 0.79 to 2.03 under
balanced deletion.

The deposited metadata identify six distinct wild-type specimens. Tissue,
library strategy and selection, paired layout, platform, and instrument are
uniform across all six. The [GEO series record for
GSE161274](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161274)
also reports a common TruSeq Stranded Total RNA and HiSeq 2000 workflow.

However, sex is absent from every deposited sample record, and no explicit
library-preparation batch field is available. Because the run accessions form
consecutive P12 and P20 blocks, an unrecorded stage-correlated batch cannot be
excluded. The robustness analysis rules out dependence on a single observed
animal; it cannot remove unmeasured confounding shared by a whole stage.

## Total-expression normalization

The source study's deposited Trinity/RSEM gene-count matrix provides an
independent total-expression check. Exact DNMT3A sequence matching against
the original source-study assembly identifies fragmented DNMT3A-like gene
components whose pooled CPM is 59.806 at P12 and 40.475 at P20, a 1.48-fold
difference. All P12 animals exceed all P20 animals (exact one-sided
permutation P = 0.05).

After dividing the exact internal-junction rate by this gene-level CPM, the
P12/P20 fold remains 2.30, with complete animal separation and P = 0.05.
Higher total DNMT3A-like expression therefore explains part, but not all, of
the internal-junction difference. The result is stable across exact-31-mer
candidate thresholds from 100 to 1000. See
`docs/DUNNART_TOTAL_EXPRESSION.md`.

## Competitive transcript-quantification sensitivity analysis

A Salmon selective-alignment analysis tested whether competitive
whole-transcript assignment supports the exact-junction result. The reference
contained all ten dunnart RefSeq DNMT3A models plus 2,093,960 non-DNMT3A
Trinity transcripts as decoys. Obvious DNMT3A Trinity representations were
identified by an explicit exact-31-mer audit and excluded from the decoy set.

The pooled internal-`066/067` estimate was 11.7962 fragments per million
input at P12 and 5.8912 at P20, a 2.00-fold difference. The internal family
represented 42.15% of assigned DNMT3A fragments at P12 and 31.62% at P20, a
1.33-fold difference. All three P12 animals exceeded all three P20 animals
for the absolute rate (exact one-sided permutation P = 0.05), but the
internal-fraction ranges overlap (P = 0.10).

However, the internal/full-length transcript-family ratio overlapped between
stages (P12 range 3.38–3.82; P20 range 1.78–4.04; P = 0.20). Competitive
quantification therefore independently supports higher internal-model signal
at P12, but the exact diagnostic junction analysis remains the stronger
evidence for increased internal-transcript signal. Neither method establishes
an internal-over-full-length isoform switch. Closely related models,
three biological replicates per stage, and the absence of capped 5′ data remain
important limitations. Full methods and outputs are documented in
`docs/DUNNART_ISOFORM_QUANTIFICATION.md`.

## Reproducible outputs

All outputs are under `results/06_isoform_evolution`:

- `marsupial_sra_inventory.tsv`
- `marsupial_sra_priority_targets.tsv`
- `dunnart_transcriptome_first_exon_evidence.tsv`
- `dunnart_transcriptome_first_exon_summary.json`
- `tammar_ont_first_exon_evidence.tsv`
- `tammar_ont_junction_read_alignments.tsv`
- `tammar_ont_first_exon_summary.json`
- `dunnart_cerebellum_ont_first_exon_evidence.tsv`
- `dunnart_cerebellum_ont_isoform_classification.tsv`
- `dunnart_cerebellum_ont_isoform_summary.json`
- `dunnart_P20_neocortex_rep1_junction_evidence.tsv`
- `dunnart_P20_neocortex_rep2_junction_evidence.tsv`
- `dunnart_P20_neocortex_rep3_junction_evidence.tsv`
- `dunnart_P20_neocortex_rep1_junction_seed_specificity.tsv`
- `dunnart_P20_neocortex_rep2_junction_seed_specificity.tsv`
- `dunnart_P20_neocortex_junction_replicates.tsv`
- `dunnart_P20_neocortex_junction_replicates_summary.json`
- `dunnart_P12_neocortex_rep2_junction_evidence.tsv`
- `dunnart_P12_neocortex_rep3_junction_evidence.tsv`
- `dunnart_P12_neocortex_rep1_junction_evidence.tsv`
- `dunnart_neocortex_stage_junctions.tsv`
- `dunnart_neocortex_stage_junctions_summary.json`
- `dunnart_neocortex_sample_metadata.tsv`
- `dunnart_neocortex_single_animal_deletion.tsv`
- `dunnart_neocortex_balanced_deletion.tsv`
- `dunnart_neocortex_stage_robustness_summary.json`
- `dunnart_source_dnmt3a_candidate_audit.tsv`
- `dunnart_source_dnmt3a_expression.tsv`
- `dunnart_source_dnmt3a_expression_summary.json`
- `dunnart_salmon_model_estimates.tsv`
- `dunnart_salmon_family_estimates.tsv`
- `dunnart_salmon_replicate_metrics.tsv`
- `dunnart_salmon_isoform_summary.json`

Rebuild with:

```bash
python scripts/inventory_marsupial_sra.py
python scripts/dunnart_transcriptome_first_exon.py
python scripts/tammar_ont_first_exon.py
python scripts/dunnart_cerebellum_ont_first_exon.py
python scripts/dunnart_cerebellum_isoform_classification.py
python scripts/dunnart_neocortex_junction_scan.py \
  --run SRR13036046 --sample P20_neocortex_rep1 \
  --sample-description "P20 neocortex replicate 1, specimen 666" \
  --output-prefix dunnart_P20_neocortex_rep1_junction
python scripts/dunnart_neocortex_junction_scan.py \
  --run SRR13036047 --sample P20_neocortex_rep2 \
  --sample-description "P20 neocortex replicate 2, specimen 746" \
  --output-prefix dunnart_P20_neocortex_rep2_junction
python scripts/summarize_dunnart_neocortex_junctions.py
python scripts/dunnart_neocortex_junction_scan.py \
  --run SRR13036043 --sample P12_neocortex_rep1 \
  --sample-description "P12 neocortex replicate 1, specimen 436" \
  --output-prefix dunnart_P12_neocortex_rep1_junction
python scripts/dunnart_neocortex_junction_scan.py \
  --run SRR13036044 --sample P12_neocortex_rep2 \
  --sample-description "P12 neocortex replicate 2, specimen 451" \
  --output-prefix dunnart_P12_neocortex_rep2_junction
python scripts/dunnart_neocortex_junction_scan.py \
  --run SRR13036045 --sample P12_neocortex_rep3 \
  --sample-description "P12 neocortex replicate 3, specimen 458" \
  --output-prefix dunnart_P12_neocortex_rep3_junction
python scripts/dunnart_neocortex_junction_scan.py \
  --run SRR13036048 --sample P20_neocortex_rep3 \
  --sample-description "P20 neocortex replicate 3, specimen 758" \
  --output-prefix dunnart_P20_neocortex_rep3_junction
python scripts/summarize_dunnart_neocortex_stages.py
python scripts/dunnart_stage_robustness.py
python scripts/dunnart_source_dnmt3a_expression.py
python scripts/prepare_dunnart_salmon_reference.py
# Build the Salmon index and quantify the six libraries as documented in
# docs/DUNNART_ISOFORM_QUANTIFICATION.md.
python scripts/summarize_dunnart_salmon_isoforms.py
python -m unittest discover -s tests
```
