# Gli3-Sufu Peak Analysis Toolkit

A comprehensive toolkit for integrating and analyzing ChIP-seq, CUT&RUN, ATAC-seq, and Split-DamID data to study gene regulatory networks, with a focus on the Gli3-Sufu pathway.

## Overview

This repository contains scripts and workflows for processing and integrating multiple types of genomic data:

- **ChIP-seq**: For identifying transcription factor binding sites
- **CUT&RUN**: For high-resolution protein-DNA interaction mapping
- **ATAC-seq**: For identifying regions of open chromatin
- **Split-DamID**: For mapping protein-DNA interactions in vivo without antibodies

The workflows automate the process from raw sequencing data processing to peak calling, differential analysis, annotation, and integration of different data types.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [ChIP-seq Analysis](#chip-seq-analysis)
- [ATAC-seq Analysis](#atac-seq-analysis)
- [CUT&RUN Analysis](#cutrun-analysis)
- [Split-DamID Analysis](#split-damid-analysis)
- [Data Integration](#data-integration)
- [License](#license)

## Installation

### Prerequisites

- Bash (>= 4.0)
- Python (>= 3.8)
- R (>= 4.0)
- Conda/Mamba

### Setting up the environment

Clone the repository:

```bash
git clone https://github.com/cfchang63/Gli3Sufu.git
cd Gli3Sufu
```

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate gli3-sufu
```

## Quick Start

The repository provides a master workflow script that can orchestrate analyses for all data types:

```bash
./scripts/gli3-sufu-master-workflow.sh \
  --output results \
  --genome mm10 \
  --chipseq-samples samples/chipseq_samples.csv \
  --atacseq-samples samples/atacseq_samples.csv \
  --cutrun-samples samples/cutrun_samples.csv \
  --splitdamid-samples samples/splitdamid_samples.txt
```

For individual analyses, see the respective sections below.

## ChIP-seq Analysis

The ChIP-seq workflow processes raw FASTQ files through alignment, peak calling, and annotation.

### Sample Preparation

Create a sample sheet CSV file with the following format:

```csv
sample,fastq_1,fastq_2,antibody,control
WT_Input,/path/to/WT_Input_R1.fastq.gz,/path/to/WT_Input_R2.fastq.gz,Input,
WT_Gli3_1,/path/to/WT_Gli3_1_R1.fastq.gz,/path/to/WT_Gli3_1_R2.fastq.gz,Gli3,WT_Input
```

### Running the Analysis

```bash
./scripts/chipseq-nfcore-workflow.sh \
  --sample-sheet samples/chipseq_samples.csv \
  --output chipseq_results \
  --genome mm10
```

### Additional Options

- `--broad-peak`: Call broad peaks instead of narrow peaks
- `--single-end`: Process single-end sequencing data
- `--integration-peaks FILE`: Integrate with external peak files

### Output

The workflow produces:
- Aligned BAM files
- Peak calls in BED format
- Peak annotations
- Genome browser tracks (BigWig format)

## ATAC-seq Analysis

The ATAC-seq workflow processes raw FASTQ files to identify open chromatin regions.

### Sample Preparation

Create a sample sheet CSV file with the following format:

```csv
sample,fastq_1,fastq_2,replicate,single_end
WT_1,/path/to/WT_1_R1.fastq.gz,/path/to/WT_1_R2.fastq.gz,1,0
WT_2,/path/to/WT_2_R1.fastq.gz,/path/to/WT_2_R2.fastq.gz,2,0
```

### Running the Analysis

```bash
./scripts/atacseq-nfcore-workflow.sh \
  --sample-sheet samples/atacseq_samples.csv \
  --output atacseq_results \
  --genome mm10
```

### Additional Options

- `--chip-peaks DIR`: Integrate with ChIP-seq peaks
- `--cutrun-peaks DIR`: Integrate with CUT&RUN peaks

### Output

The workflow produces:
- Aligned BAM files
- Peak calls in BED format
- Differential accessibility analysis
- Genome browser tracks (BigWig format)
- Integration with other data types

## CUT&RUN Analysis

CUT&RUN provides high-resolution mapping of protein-DNA interactions with lower background than ChIP-seq.

### Sample Preparation

Create a sample sheet CSV file with the following format (similar to ChIP-seq):

```csv
sample,fastq_1,fastq_2,antibody,control
WT_Input,/path/to/WT_Input_R1.fastq.gz,/path/to/WT_Input_R2.fastq.gz,Input,
WT_Gli3_1,/path/to/WT_Gli3_1_R1.fastq.gz,/path/to/WT_Gli3_1_R2.fastq.gz,Gli3,WT_Input
```

### Running the Analysis

```bash
./scripts/cutandrun-peak-intersect-annotate.sh \
  --output cutrun_results \
  --genome mm10 \
  --run-nf-core \
  --sample-sheet samples/cutrun_samples.csv \
  --workflow cutandrun
```

### Output

The workflow produces:
- Aligned BAM files
- Peak calls using SEACR (specialized for CUT&RUN)
- Genome browser tracks (BigWig format)
- Peak annotations

## Split-DamID Analysis

Split-DamID maps protein-DNA interactions in vivo without antibodies or fixation.

### Sample Preparation

Create a text file listing your samples:

```
DAM-1
DAM-2
DAM-3
DAM-Dox-1
DAM-Dox-2
DAM-Dox-3
Gli3-Hand2-1
Gli3-Hand2-2
Gli3-Hand2-3
Gli3-Hand2-Dox-1
Gli3-Hand2-Dox-2
Gli3-Hand2-Dox-3
```

### Preparing GATC Sites and Blacklist

Before running Split-DamID analysis, you need to prepare:

1. GATC sites file:
```bash
python scripts/create-gatc-sites.py /path/to/genome.fa /path/to/gatc_sites.bed
```

2. Blacklist + mtDNA file:
```bash
./scripts/create-blacklist-mtdna.sh --genome mm10 --output blacklist_mtDNA.bed
```

### Running the Analysis

```bash
./scripts/splitdamid-workflow.sh \
  --samples samples/splitdamid_samples.txt \
  --output splitdamid_results \
  --genome mm10 \
  --fasta /path/to/mm10.fa \
  --bowtie2 /path/to/mm10_bt2index \
  --gatc /path/to/gatc_sites.bed \
  --blacklist /path/to/blacklist_mtDNA.bed \
  --extension-script scripts/gatc-extension-script.py
```

### Differential Analysis

After processing, run differential binding analysis:

```bash
./scripts/splitdamid-binning-workflow.sh \
  --input splitdamid_results \
  --output splitdamid_diffbind
```

### Output

The workflow produces:
- BAM files (aligned, extended to GATC sites, filtered)
- Peak calls from MACS3
- Coverage tracks (BigWig, bedGraph)
- Binned genomic coverage analysis
- Log2 fold change calculations between conditions

## Data Integration

The toolkit provides scripts for integrating various data types:

```bash
./scripts/peak-intersection-workflow.sh \
  --output integration_results \
  --genome mm10 \
  --chipseq-peaks chipseq_results/peaks \
  --atacseq-peaks atacseq_results/peaks \
  --cutrun-peaks cutrun_results/peaks \
  --splitdamid-peaks splitdamid_results/peaks
```

For motif analysis of intersections:

```bash
./scripts/motif-analysis-homer.sh \
  --output motif_results \
  --genome mm10 \
  --intersections integration_results/intersections
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
