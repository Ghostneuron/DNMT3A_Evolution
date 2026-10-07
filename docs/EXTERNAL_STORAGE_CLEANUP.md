# External storage cleanup log

## 2026-07-26: checksum-failed FASTQ copies

To preserve space for the active full mCH analysis and another external-drive
job, the following checksum-failed, scientifically unusable copies were
deleted after checksum-valid replacements had been obtained:

- `fastq/SRR16892996_1.fastq.gz.corrupt-md5-3d73c06a`
- `fastq/SRR16892996_2.fastq.gz.corrupt-md5-0822674e`
- `quarantine/failed_md5_20260725/SRR16892995_2.fastq.gz`
- `quarantine/failed_md5_20260725/SRR16892997_1.fastq.gz`
- `quarantine/failed_md5_20260725/SRR16892997_2.fastq.gz`
- `quarantine/failed_md5_20260725/SRR16893002_2.fastq.gz`
- `quarantine/failed_md5_20260725_attempt2/SRR16892997_1.fastq.gz`

Associated `.aria2` state files were also removed. These files were
regenerable from ENA but were not valid scientific inputs. All 16
checksum-valid FASTQs were retained.

## 2026-08-02: validated downstream-successor cleanup

The storage-aware cleanup tool validated the retained deduplicated BAMs with
`samtools quickcheck` before deleting any predecessor. For `SRR16892995`, the
completed coverage, M-bias, splitting report, context table, autosomal-bin
table, and gene-body table were also verified before removing the two trimmed
FASTQs, three redundant Bismark context files, and redundant bedGraph. This
reclaimed 42.724 GiB. The original checksum-valid raw FASTQs, deduplicated BAM,
CX coverage file, compact summaries, and QC reports remain available.

For `SRR16892996`, the validated deduplicated BAM superseded the alignment BAM;
removing that alignment BAM reclaimed 20.847 GiB. Its other intermediates are
being retained until context summarization completes.

An inventory on 2026-08-02 found that `SRR16892997_2.fastq.gz` and
`SRR16893002_2.fastq.gz` were no longer present in the active raw directory,
despite the earlier complete-set validation. Checksum-enforced, single-stream
restoration from ENA was started before either affected library was processed.

## 2026-08-06: WT2 and KO1 successor cleanup

After `SRR16892996` context, autosomal-bin, and gene-body summaries were
validated, its trimmed FASTQs, redundant three-context files, and redundant
bedGraph were removed, reclaiming 34.840 GiB. Its raw FASTQs, deduplicated BAM,
CX coverage, M-bias data, reports, and compact summaries remain available.

The `SRR16892997` deduplicated BAM passed `samtools quickcheck`, allowing its
superseded alignment BAM to be removed and reclaiming 21.037 GiB. Its remaining
intermediates are retained until the running context summarization validates.

## 2026-08-11: KO2 summarized-successor cleanup

After `SRR16892998` context, autosomal-bin, and gene-body summaries were
validated, its superseded alignment BAM, trimmed FASTQs, redundant
three-context files, and redundant bedGraph were removed. This reclaimed
59.765 GiB while retaining the raw FASTQs, deduplicated BAM, CX coverage,
M-bias data, QC reports, and compact summaries.

## 2026-08-18: final four-library successor cleanup

After each library's context, autosomal-bin, and gene-body summaries and the
required retained successors were validated, redundant alignment and
extraction intermediates were removed. Reclaimed space was:

- `SRR16892999`: 63.290 GiB
- `SRR16893000`: 52.842 GiB
- `SRR16893001`: 56.217 GiB
- `SRR16893002`: 52.048 GiB

The cleanup retained the scientific inputs and validated downstream products,
including raw FASTQs, deduplicated BAMs, CX coverage, M-bias and QC reports,
and compact context summaries. The completed eight-library comparison and
final replicate-level report were then generated. External-drive free space
was approximately 592 GiB after cleanup.
