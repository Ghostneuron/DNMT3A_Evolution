# DNMT3A1 versus DNMT3A2 neuronal mCH reanalysis

## Completion status (2026-08-18)

The full eight-library workflow is complete, with two replicates for every
genotype. The final audit and interpretation are in `docs/MCH_FINAL_RESULTS.md`.
Machine-readable replicate and genotype summaries are in
`results/bismark/final_summary/` on the external analysis drive.

The primary result is that `Dnmt3a1_KO` reduces lambda-corrected autosomal mCA
by approximately 99.4% relative to WT, whereas `Dnmt3a2_KO` retains bulk mCA.
The `Dnmt3a1_delta_N` samples retain approximately 40% of WT corrected mCA,
supporting a substantial contribution from the DNMT3A1 N-terminal region.
These are replicated descriptive effects (n=2 per genotype), not strong
population-level statistical estimates.

## Question

Does selective loss of DNMT3A1 or DNMT3A2 produce distinguishable effects on
postnatal neuronal non-CpG methylation, especially mCA?

The eight P21 NeuN-positive cortical-neuron EM-seq libraries in GSE164265
provide two biological replicates each for wild type, `Dnmt3a1` knockout,
`Dnmt3a2` knockout, and a `Dnmt3a1` N-terminal deletion. The deposited
processed tracks are CpG-only, so the raw paired-end reads must be recalled in
all cytosine contexts.

## Why this is a new analysis of an existing experiment

The source study's public WGBS workflow mapped reads with BSMAP and then called
methylation with `methratio.py -x CG`; that option restricted the reported
calls to CpG. Its downstream DMR and profile analyses therefore addressed CpG
methylation, not neuronal mCH. The raw reads permit a distinct all-context
reanalysis. The source workflow is archived at
<https://github.com/DapengHao/DNMT3A1N_terminus/blob/main/WGBS%20analysis.md>.

## Storage and provenance

Large files and analysis products are isolated on:

`$DNMT3A_MCH_ROOT` (or the repository-relative default `external_data/mCH`)

The local manifest `config/mch_samples.tsv` records the GEO/SRA mapping, exact
ENA FASTQ URLs, byte counts, and MD5 checksums. `scripts/download_mch_fastq.py`
resumes downloads and accepts a file only after both size and checksum match.
The 16 paired FASTQ files total 147.30 GB compressed.
The reference is UCSC mm9, matching the assembly reported by the source study,
with the official NEB lambda and pUC19 control sequences appended.
The composite retains UCSC random contigs to reduce spurious placement of
repeat-derived reads; the source repository removed `*_random.fa`. Primary
statistics are restricted to chr1–chr19, and any material difference from the
deposited CpG tracks will trigger a primary-chromosome-only sensitivity run.

## Computational workflow

1. Build Bismark/Bowtie2 bisulfite indexes from mm9 plus the EM-seq controls.
2. Detect and trim paired-end adapters without quality- or poly-G-based read
   removal, which could distort low-complexity converted reads.
3. Pilot one wild-type replicate on a fixed number of read pairs.
4. Confirm mapping mode, duplicate rate, context counts, and per-cycle M-bias.
5. Choose end-trimming from M-bias profiles rather than assuming zero trimming.
6. Align and deduplicate all eight libraries.
7. Extract CpG, CHG, and CHH calls with overlapping mate calls counted once.
8. Separate mCA, mCT, and mCC from the sequence context in the mm9 reference.
9. Compare genotype effects using biological replicates and report uncertainty.

The pipeline is `scripts/run_mch_bismark.py`. Its extraction stage uses
`--CX --comprehensive --no_overlap` and deliberately does not merge CHG with
CHH. It also deliberately avoids Bismark's non-conversion filter, because
clusters of methylated CH are expected biological signal in neurons.
`scripts/summarize_mch_contexts.py` reorients calls on both DNA strands,
classifies CA, CT, CC, CG, CHG, and CHH, and produces autosomal 100-kb and
optional gene-body summaries. `scripts/compare_mch_genotypes.py` then reports replicate-level
values, corrects non-CpG fractions with the matched lambda control when it is
present, and reports each mutant-versus-WT contrast without treating genomic
bins as independent biological replicates.
The gene-body reference contains 21,279 longest protein-coding RefSeq
transcripts prepared from the UCSC mm9 `refGene` table.

## Required quality controls

- FASTQ size and MD5 validation.
- Mapping rate and bisulfite strand balance.
- PCR duplicate fraction.
- Per-cycle M-bias, separately for R1 and R2.
- Coverage and call counts for CpG, CA, CT, and CC.
- Conversion-error estimate from unmethylated lambda and protection efficiency
  from CpG-methylated pUC19, if the expected spike-in reads are present.
- Concordance between biological replicates.
- CpG means compared with the deposited CpG tracks as a pipeline sanity check.

## Checksum-validated WT pilot

WT replicate 2 (`SRR16892996`) passed independent size and ENA MD5 validation.
A fastp pass over its first 100,000 read pairs read all 200,000 mates
successfully. Before trimming, 89.62% of bases were Q30 and mean read length was
151 bp. Adapter evidence was detected in 22,912 reads (11.46%); 675,010 bases
were removed, leaving mean R1/R2 lengths of 147/148 bp.

Bismark assigned unique best alignments to 59,668/100,000 pairs. Deduplication
retained 55,959 alignments (93.78%; 6.22% duplicates). The untrimmed M-bias
diagnostic showed elevated non-CpG calls at the R1 5-prime end: pooled CHH was
2.43% at bases 1--5 and 2.08% at bases 6--10, versus 1.52% at bases 11--140.
R1 bases 141--151 were modestly elevated (1.72%). R2 non-CpG rates were stable,
although its first five bases showed a CpG-specific distortion. The fixed
extraction preset is therefore R1 10 bp at both ends and R2 5 bp at the
5-prime end.

After this trimming, the pilot contained 2,862,482 classified cytosine rows
and no unknown contexts. Autosomal fractions were CpG 66.59%, CA 2.41%, CC
0.89%, and CT 1.01%. The deposited CpG-only value for the same replicate is
67.67%, providing a close pipeline sanity check at pilot depth. Unmethylated
lambda had 0/532 methylated CH calls; pUC19 had only 21 CpG calls in the small
pilot, so protection efficiency will be judged from full libraries.

## Primary and secondary tests

The primary endpoint is the genotype effect on autosomal mCA. The
`Dnmt3a1_delta_N` contrast tests whether the isoform-specific N-terminal tail,
rather than only total DNMT3A1 abundance, is required for the methylation
pattern. Aggregate mCH is secondary because neuronal mCA is the dominant
biologically interpretable component. Further analyses will include gene-body
mCA, long genes, MeCP2-sensitive genes, and genomic regions associated with
H3K36me2 and H2AK119ub.

The three pre-specified primary contrasts are `Dnmt3a1_KO - WT`,
`Dnmt3a2_KO - WT`, and `Dnmt3a1_delta_N - WT`. Effect sizes and agreement
between the two biological replicates are primary; with only two animals per
group, nominal site/bin-level p-values will not be presented as strong
biological replication. Regional analyses will use chromosome-block
resampling or an explicitly overdispersed methylation model and will remain
secondary to the replicate-level global result.

This experiment can establish an isoform-specific association with neuronal
mCH deposition. It cannot by itself show that DNMT3A2 caused the evolutionary
expansion of brain size or intelligence.

## Reproducible entry points

```bash
python scripts/download_mch_fastq.py
python scripts/prepare_mch_reference.py
python scripts/run_mch_bismark.py prepare --threads 2
python scripts/run_mch_bismark.py all --sample SRR16892996 \
  --upto 100000 --threads 2 --ignore 10 --ignore-r2 5 \
  --ignore-3prime 10
python scripts/summarize_mch_contexts.py SAMPLE.bismark.cov.gz \
  --output-dir SAMPLE_CONTEXT_DIRECTORY
python scripts/compare_mch_genotypes.py
```

The `--upto` command is a pilot only. Full runs omit it and retain the
M-bias-derived extraction preset shown above. Per-sample M-bias profiles will
still be audited; any genotype-correlated deviation will trigger a uniform
sensitivity trimming rule rather than genotype-specific trimming.
