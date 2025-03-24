# ChIP-seq Analysis Workflow

This document provides detailed information about the ChIP-seq analysis workflow implemented in this repository.

## Overview

ChIP-seq (Chromatin Immunoprecipitation followed by sequencing) is a method used to analyze protein interactions with DNA. This workflow processes raw ChIP-seq data through alignment, peak calling, and annotation to identify genome-wide binding sites for proteins of interest.

## Workflow Steps

The ChIP-seq workflow consists of the following steps:

1. **Quality Control**: Raw reads are checked for quality using FastQC
2. **Alignment**: Reads are aligned to the reference genome using Bowtie2
3. **Post-alignment Processing**: 
   - BAM filtering
   - Duplicate removal
   - Blacklist filtering
4. **Peak Calling**: MACS2 is used to identify enriched regions
5. **Peak Annotation**: HOMER is used to annotate peaks with genomic features
6. **Visualization**: BigWig files are generated for visualization

## Required Input

### Sample Sheet Format

The workflow requires a CSV sample sheet with the following format:

```csv
sample,fastq_1,fastq_2,antibody,control
WT_Input,/path/to/WT_Input_R1.fastq.gz,/path/to/WT_Input_R2.fastq.gz,Input,
WT_Gli3_1,/path/to/WT_Gli3_1_R1.fastq.gz,/path/to/WT_Gli3_1_R2.fastq.gz,Gli3,WT_Input
```

Fields:
- `sample`: Unique sample identifier
- `fastq_1`: Path to R1 FASTQ file
- `fastq_2`: Path to R2 FASTQ file (leave empty for single-end data)
- `antibody`: Antibody used for IP (e.g., Gli3, SATB2)
- `control`: Control sample identifier (usually Input)

### Reference Files

- Reference genome FASTA file
- Genome index for Bowtie2
- Blacklist regions in BED format
- Chromosome sizes file

## Running the Workflow

### Basic Usage

```bash
./scripts/chipseq-nfcore-workflow.sh \
  --sample-sheet samples/chipseq_samples.csv \
  --output chipseq_results \
  --genome mm10
```

### Full Options

```bash
./scripts/chipseq-nfcore-workflow.sh \
  --sample-sheet samples/chipseq_samples.csv \
  --output chipseq_results \
  --genome mm10 \
  --genome-dir /path/to/genome_dir \
  --fasta /path/to/genome.fa \
  --gtf /path/to/annotation.gtf \
  --blacklist /path/to/blacklist.bed \
  --bowtie2-index /path/to/bowtie2_index \
  --macs-gsize 2.7e9 \
  --threads 16 \
  --single-end \
  --read-length 75 \
  --broad-peak \
  --integration-peaks /path/to/external_peaks.bed
```

## Output Files

The workflow generates the following output directories:

- `chipseq_results/`: Main output directory
  - `nf-core-chipseq/`: nf-core pipeline output
    - `fastqc/`: Quality control reports
    - `bowtie2/`: Alignment results
      - `mergedLibrary/`: Merged replicates
        - `macs2/`: Peak calling results
          - `narrowPeak/` or `broadPeak/`: Peak files
          - `consensus/`: Consensus peaks across replicates
  - `bam/`: Processed BAM files
  - `bigwig/`: BigWig files for visualization
  - `peaks/`: Peak files in BED format
  - `annotation/`: Annotated peak files
  - `integration/`: Integration with other data types (if specified)

### Key Output Files

- `*_peaks.narrowPeak` or `*_peaks.broadPeak`: Called peaks
- `*.consensus_peaks.bed`: Consensus peaks across replicates
- `*_annotation.txt`: HOMER annotation of peaks
- `*.bw`: BigWig files for visualization

## Integration with Other Data Types

The ChIP-seq workflow can integrate with other data types by using the `--integration-peaks` option:

```bash
./scripts/chipseq-nfcore-workflow.sh \
  --sample-sheet samples/chipseq_samples.csv \
  --output chipseq_results \
  --genome mm10 \
  --integration-peaks atacseq_results/consensus_peaks.bed
```

This will identify overlapping regions between ChIP-seq peaks and the provided peak file.

## Downstream Analysis

After running the ChIP-seq workflow, you can:

1. **Annotate Peaks**:
```bash
Rscript scripts/chipseq-annotation.r \
  --input chipseq_results/peaks/Sample_peaks.narrowPeak \
  --output chipseq_results/annotation/Sample_annotation.txt \
  --genome mm10 \
  --tss-region -3000,3000 \
  --go-analysis
```

2. **Integrate with Other Data**:
```bash
./scripts/peak-intersection-workflow.sh \
  --output integration_results \
  --genome mm10 \
  --chipseq-peaks chipseq_results/peaks
```

3. **Perform Motif Analysis**:
```bash
./scripts/motif-analysis-homer.sh \
  --output motif_results \
  --genome mm10 \
  --intersections integration_results/intersections
```

## Troubleshooting

### Common Issues

1. **Missing input files**: Ensure all FASTQ files exist at the specified paths
2. **Index compatibility**: Make sure the genome index matches the reference genome
3. **Peak calling failure**: Check if control samples are properly specified

### Logs and Reports

- Check log files in the output directory for errors
- Examine MultiQC reports for quality metrics
- Review nf-core pipeline execution reports

## References

- nf-core ChIP-seq pipeline: [https://nf-co.re/chipseq](https://nf-co.re/chipseq)
- HOMER: [http://homer.ucsd.edu/homer/](http://homer.ucsd.edu/homer/)
- MACS2: [https://github.com/macs3-project/MACS](https://github.com/macs3-project/MACS)
