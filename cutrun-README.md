# CUT&RUN Analysis Workflow

This document provides detailed information about the CUT&RUN analysis workflow implemented in this repository.

## Overview

CUT&RUN (Cleavage Under Targets and Release Using Nuclease) is a chromatin profiling method that provides high-resolution mapping of protein-DNA interactions with lower background than ChIP-seq. This workflow processes raw CUT&RUN data through alignment, peak calling, and annotation to identify genome-wide binding sites for proteins of interest.

## Workflow Steps

The CUT&RUN workflow consists of the following steps:

1. **Quality Control**: Raw reads are checked for quality using FastQC
2. **Adapter Trimming**: Adapters are removed with Trimmomatic
3. **Alignment**: Reads are aligned to the reference genome using Bowtie2
4. **Post-alignment Processing**:
   - BAM filtering
   - Blacklist filtering
   - Spike-in normalization (if available)
5. **Peak Calling**: SEACR (specialized for CUT&RUN) and/or MACS2 are used to identify enriched regions
6. **Peak Annotation**: HOMER is used to annotate peaks with genomic features
7. **Visualization**: BigWig files are generated for visualization

## Required Input

### Sample Sheet Format

The workflow requires a CSV sample sheet with the following format (similar to ChIP-seq):

```csv
sample,fastq_1,fastq_2,antibody,control
WT_IgG,/path/to/WT_IgG_R1.fastq.gz,/path/to/WT_IgG_R2.fastq.gz,IgG,
WT_Gli3_1,/path/to/WT_Gli3_1_R1.fastq.gz,/path/to/WT_Gli3_1_R2.fastq.gz,Gli3,WT_IgG
```

Fields:
- `sample`: Unique sample identifier
- `fastq_1`: Path to R1 FASTQ file
- `fastq_2`: Path to R2 FASTQ file
- `antibody`: Antibody used for IP (e.g., Gli3, SATB2)
- `control`: Control sample identifier (usually IgG)

### Reference Files

- Reference genome FASTA file
- Genome index for Bowtie2
- Blacklist regions in BED format
- Chromosome sizes file
- Spike-in genome (optional, e.g., E. coli)

## Running the Workflow

### Basic Usage

```bash
./scripts/cutandrun-peak-intersect-annotate.sh \
  --output cutrun_results \
  --genome mm10 \
  --run-nf-core \
  --sample-sheet samples/cutrun_samples.csv \
  --workflow cutandrun
```

### Full Options

```bash
./scripts/cutandrun-peak-intersect-annotate.sh \
  --output cutrun_results \
  --genome mm10 \
  --gtf /path/to/annotation.gtf \
  --fasta /path/to/genome.fa \
  --blacklist /path/to/blacklist.bed \
  --bowtie2 /path/to/bowtie2_index \
  --homer /path/to/homer/annotatePeaks.pl \
  --run-nf-core \
  --sample-sheet samples/cutrun_samples.csv \
  --workflow cutandrun
```

## Output Files

The workflow generates the following output directories:

- `cutrun_results/`: Main output directory
  - `nf-core-cutandrun/`: nf-core pipeline output
    - `fastqc/`: Quality control reports
    - `trimgalore/`: Trimmed reads
    - `bowtie2/`: Alignment results
    - `seacr/`: SEACR peak calling results
    - `macs2/`: MACS2 peak calling results (if used)
    - `bigwig/`: BigWig coverage files
  - `bed_files/`: Processed BED files
  - `intersections/`: Peak intersection results
  - `annotations/`: Annotated peak files

### Key Output Files

- `*_peaks.stringent.bed`: SEACR stringent peaks
- `*_peaks.relaxed.bed`: SEACR relaxed peaks
- `*.consensus_peaks.bed`: Consensus peaks across replicates
- `*_annotation.txt`: HOMER annotation of peaks
- `*.bw`: BigWig files for visualization
- `*.overlap.bed`: Intersection results with other data types

## CUT&RUN-Specific Considerations

CUT&RUN has several specific considerations compared to ChIP-seq:

### Specialized Peak Calling with SEACR

SEACR (Sparse Enrichment Analysis for CUT&RUN) is specifically designed for CUT&RUN data:

```bash
# SEACR peak calling
seacr_path=/path/to/SEACR/SEACR_1.3.sh
bash ${seacr_path} ${target_bedgraph} ${control_bedgraph} norm stringent ${output_prefix}
```

The workflow supports both stringent and relaxed thresholds, with stringent being more selective.

### Spike-in Normalization

If spike-in DNA (e.g., E. coli) is available:

```bash
# Calculate spike-in normalization factor
spike_in_ratio=$(samtools view -c -F 260 ${bam_file} ${spike_in_chromosome} | awk -v total=${total_reads} '{print total/$1}')
```

### Fragment Size Distribution

CUT&RUN typically produces smaller fragments than ChIP-seq:

- Fragment sizes of ~180 bp indicate nucleosome-protected regions
- Fragment sizes of <120 bp indicate transcription factor binding sites

## Integration with Other Data Types

The CUT&RUN workflow can integrate with ChIP-seq, ATAC-seq, and Split-DamID data:

```bash
./scripts/peak-intersection-workflow.sh \
  --output integration_results \
  --genome mm10 \
  --cutrun-peaks cutrun_results/seacr \
  --chipseq-peaks chipseq_results/peaks \
  --atacseq-peaks atacseq_results/peaks
```

This identifies overlapping regions between different data types to discover multi-evidence binding sites.

## Downstream Analysis

After running the CUT&RUN workflow, you can:

1. **Annotate Peaks**:
```bash
./scripts/cutandrun-peak-intersect-annotate.sh \
  --output cutrun_annotation \
  --genome mm10 \
  --homer /path/to/homer/annotatePeaks.pl \
  --bed-files cutrun_results/bed_files
```

2. **Perform Motif Analysis**:
```bash
./scripts/motif-analysis-homer.sh \
  --output motif_results \
  --genome mm10 \
  --intersections cutrun_results/intersections
```

3. **Generate Heatmaps and Aggregate Plots**:
```bash
computeMatrix reference-point -S sample.bw -R peaks.bed -a 2000 -b 2000 -o matrix.gz
plotHeatmap -m matrix.gz -o heatmap.png
plotProfile -m matrix.gz -o profile.png
```

## Troubleshooting

### Common Issues

1. **High background signal**: Check IgG controls and SEACR parameters
2. **Fragment size distribution issues**: Check if MNase or other enzymatic digestion worked properly
3. **Low signal-to-noise ratio**: Adjust peak calling parameters or increase sequencing depth

### Quality Metrics

Key quality metrics for CUT&RUN:
- FRiP (Fraction of Reads in Peaks)
- Fragment size distribution
- Enrichment over IgG control
- Specificity (signal at expected targets)

## References

- nf-core CUT&RUN pipeline: [https://nf-co.re/cutandrun](https://nf-co.re/cutandrun)
- SEACR: [https://github.com/FredHutch/SEACR](https://github.com/FredHutch/SEACR)
- Original CUT&RUN protocol: [https://doi.org/10.1038/nprot.2018.015](https://doi.org/10.1038/nprot.2018.015)
- Improved CUT&RUN protocol: [https://doi.org/10.1038/s41467-019-09982-5](https://doi.org/10.1038/s41467-019-09982-5)