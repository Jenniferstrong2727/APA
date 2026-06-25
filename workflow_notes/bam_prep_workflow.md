# BAM Preprocessing Workflow

## Overview

The in-house single-nucleus RNA-seq data were initially provided as pooled BAM files, where each sequencing run contained reads from multiple donor-region samples. Because downstream APA analysis requires donor-specific BAM files, the pooled BAMs were split into individual donor-region BAMs using donor barcode assignments derived from the nuclei metadata.

This document describes the development of the preprocessing workflow, from an initial pilot on a single sequencing run through scaling the workflow to all runs using parallel execution on the Minerva HPC cluster.

---

# Workflow Summary

The preprocessing workflow consists of six major steps:

1. Identify donor-region samples contained within each pooled sequencing run.
2. Generate donor-specific barcode whitelist files from the nuclei metadata.
3. Normalize barcode formatting between metadata and BAM files.
4. Split pooled BAMs into donor-specific BAMs using `samtools`.
5. Index and quality-control each BAM.
6. Scale the workflow across all sequencing runs using parallel LSF jobs.

---

# Pilot Workflow (GEX-10-1)

Workflow development began with sequencing run **GEX-10-1**, which contains three donor-region samples:

* NBB20_071_AMY
* NBB19_071_AMY
* NBB19_074_AMY

The pilot workflow was used to verify:

* the correct BAM barcode tag (`CB`)
* barcode formatting
* whitelist generation
* BAM splitting
* BAM indexing
* downstream quality control

Successful completion of the pilot provided confidence before processing all remaining sequencing runs.

---

# Barcode Whitelist Generation

Barcode whitelist files are generated from the nuclei metadata after donor assignments have been resolved.

The script:

`make_nuclei_whitelists.R`

performs the following tasks:

* reads the nuclei metadata
* reads the master sample sheet
* converts metadata barcodes from `.1` to `-1`
* groups barcodes by donor-region sample
* writes one whitelist file per donor-region sample

These whitelist files are subsequently used during BAM extraction.

---

# BAM Splitting

Donor-specific BAMs are generated using:

```bash
samtools view -b -D CB:<barcode_whitelist>
```

Each output BAM is then:

* indexed with `samtools index`
* validated using `samtools quickcheck`
* checked by confirming successful read extraction

---

# Pipeline Scripts

The BAM preprocessing workflow is implemented using four primary scripts.

### make_nuclei_whitelists.R

Creates donor-specific barcode whitelist files from the nuclei metadata.

### run_split_GEX10_1.sh

Pilot script used to validate the workflow on a single sequencing run.

### run_split_one_run.sh

Generalized workflow for processing a single sequencing run using the manifest file.

Designed for submission as one LSF job per sequencing run.

### run_split_all_nuclei.sh

Processes every run listed in the manifest sequentially.

Primarily retained for testing and debugging.

---

# Executing the Scripts

Before execution, shell scripts were made executable:

```bash
chmod +x run_split_GEX10_1.sh
chmod +x run_split_one_run.sh
chmod +x run_split_all_nuclei.sh
```

The pilot workflow was executed first to verify correctness before scaling to the full dataset.

---

# HPC Job Submission

Jobs were submitted to the Minerva LSF scheduler using `bsub`.

Each submission specified:

* project allocation
* queue
* CPU allocation
* memory requirements
* wall time
* output log
* error log

Representative submission:

```bash
bsub \
-P acc_als-omics \
-q premium \
-n 1 \
-W 24:00 \
-R "span[hosts=1] rusage[mem=32000]" \
-o <stdout_log> \
-e <stderr_log> \
run_split_one_run.sh <RUN_NAME>
```

Standard output and error logs were written for every sequencing run to facilitate monitoring and troubleshooting.

---

# Parallel Processing Strategy

Following successful validation of the pilot workflow, preprocessing was scaled across all remaining sequencing runs.

Rather than processing runs sequentially, a reusable script (`run_split_one_run.sh`) was developed to process a single sequencing run.

A shell loop then submitted one independent LSF job for each sequencing run using a list of run identifiers.

This strategy allowed multiple runs to execute simultaneously across the cluster, substantially reducing total processing time while preserving independent logging and restart capability for each run.

If individual runs failed due to transient cluster or filesystem issues, only those runs were resubmitted without repeating successful jobs.

---

# Quality Control

Each donor-specific BAM underwent the following validation:

* successful BAM generation
* BAM indexing
* `samtools quickcheck`
* read count inspection
* log file review

The workflow was considered complete only after all donor-region BAMs were successfully generated and passed quality control.

---

# Repository Structure

```
scripts/
└── bam_prep/
    ├── make_nuclei_whitelists.R
    ├── run_split_GEX10_1.sh
    ├── run_split_one_run.sh
    └── run_split_all_nuclei.sh

workflow_notes/
└── bam_prep_workflow.md
```

The resulting donor-specific BAM files serve as the input for the downstream APA analysis pipeline.
