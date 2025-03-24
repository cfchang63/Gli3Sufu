# ATAC-seq Analysis Workflow

This document provides detailed information about the ATAC-seq analysis workflow implemented in this repository.

## Overview

ATAC-seq (Assay for Transposase-Accessible Chromatin with high-throughput sequencing) is a method used to assess genome-wide chromatin accessibility. This workflow processes raw ATAC-seq data through alignment, peak calling, and differential analysis to identify open chromatin regions in the genome.

## Workflow Steps

The ATAC-seq workflow consists of the following steps:

1. **Quality Control**: Raw reads are checked for quality using FastQC
2. **Adapter Trimming**: Adapters are removed with Trimmomatic
3. **Alignment**: Reads are aligned to the reference genome using Bowtie2
4. **Post-alignment Processing**:
   - Mitochondrial read removal
   - BAM filtering
   - Duplicate removal
   - Blacklist filtering
   - Read shifting (Tn5 transposase adjustment)
5. **Peak Calling**: MACS3 is used to identify open chromatin regions
6. **Differential Analysis**: Comparison between sample groups
7. **Peak Annotation**: Annotation of peaks with genomic features
8. **Visualization**: Generation of coverage tracks and enrichment profiles

## Required Input

### Sample Sheet Format

The workflow requires a CSV sample sheet with the following format:

```csv
sample,fastq_1,fastq_2,replicate,single_end
WT_1,/path/to/WT_1_R1.fastq.gz,/path/to/WT_1_R2.fastq.gz,1,0
WT_2,/path/to/WT_2_R1.fastq.gz,/path/to/WT_2_R2.fastq.gz,2,0
KO_1,/path/to/KO_1_R1.fastq.gz,/path/to/KO_1_R2.fastq.gz,1,0
KO_2,/path/to/KO_2_R1.fastq.gz,/path/to/KO_2_R2.fastq.gz,2,0
```

Fields:
- `sample`: Unique sample identifier
- `fastq_1`: Path to R1 FASTQ file
- `fastq_2`: Path to R2 FASTQ file (leave empty for single-end data)
- `replicate`: Replicate number
- `single_end`: Boolean flag for single-end data (0=paired, 1=single)

### Reference Files

- Reference genome FASTA file
- Genome index for Bowtie2
- Blacklist regions in BED format
- GTF gene annotation file
- Chromosome sizes file

## Running the Workflow

### Basic Usage

```bash
./scripts/atacseq-nfcore-workflow.sh \
  --sample-sheet samples/atacseq_samples.csv \
  --output atacseq_results \
  --genome mm10
```

### Full Options

```bash
./scripts/atacseq-nfcore-workflow.sh \
  --sample-sheet samples/atacseq_samples.csv \
  --output atacseq_results \
  --genome mm10 \
  --gtf /path/to/annotation.gtf \
  --fasta /path/to/genome.fa \
  --blacklist /path/to/blacklist.bed \
  --bowtie2 /path/to/bowtie2_index \
  --threads 16 \
  --chip-peaks /path/to/chipseq/peaks \
  --cutrun-peaks /path/to/cutrun/peaks
```

## Output Files

The workflow generates the following output directories:

- `atacseq_results/`: Main output directory
  - `nf-core-atacseq/`: nf-core pipeline output
    - `fastqc/`: Quality control reports
    - `trimgalore/`: Trimmed reads
    - `bowtie2/`: Alignment results
    - `macs/`: Peak calling results
      - `consensus/`: Consensus peaks across replicates
    - `bigwig/`: BigWig coverage files
    - `differential_accessibility/`: Differential analysis results
  - `integrations/`: Integration with other data types
  - `annotations/`: Annotated peak files

### Key Output Files

- `*_peaks.narrowPeak`: Called peaks
- `consensus_peaks.mLb.clN.bed`: Consensus peaks across replicates
- `*_annotation.txt`: Annotation of peaks
- `*.bw`: BigWig files for visualization
- `differential_peaks.txt`: Differential accessibility analysis results

## Specialized ATAC-seq Processing

The ATAC-seq workflow includes several specialized processing steps:

### Mitochondrial Read Removal

Mitochondrial reads are common contaminants in ATAC-seq data and are removed:

```bash
samtools view -h "${BAM_FILE}" | grep -v chrM | samtools sort -O bam -o "${FILTERED_BAM}"
```

### Tn5 Transposase Shift

ATAC-seq reads are shifted to account for Tn5 transposase binding:

```bash
alignmentSieve --numberOfProcessors max --ATACshift \
  --blackListFileName "${BLACKLIST}" \
  --bam "${INPUT_BAM}" -o "${SHIFTED_BAM}"
```

### Fragment Size Distribution

ATAC-seq fragment size distribution is analyzed to ensure proper nucleosome-free and nucleosome-bound patterns.

## Differential Accessibility Analysis

The workflow includes differential accessibility analysis using DESeq2:

```bash
Rscript scripts/atac-seq-diff-abundance.r \
  --counts counts_matrix.txt \
  --metadata metadata.txt \
  --output results_dir \
  --fdr 0.05 \
  --log2fc 1.0
```

This generates:
- MA plots
- Volcano plots
- PCA plots for sample clustering
- Lists of differentially accessible regions
- GO and KEGG pathway enrichment

## Integration with Other Data Types

The ATAC-seq workflow can integrate with ChIP-seq and CUT&RUN data:

```bash
./scripts/atacseq-nfcore-workflow.sh \
  --sample-sheet samples/atacseq_samples.csv \
  --output atacseq_results \
  --genome mm10 \
  --chip-peaks /path/to/chipseq/peaks \
  --cutrun-peaks /path/to/cutrun/peaks
```

This identifies overlapping regions between ATAC-seq open chromatin and transcription factor binding sites.

## Peak Annotation

ATAC-seq peaks are annotated with genomic features:

```bash
python scripts/atac-seq-peak-annotation.py \
  --input peaks.bed \
  --genome mm10 \
  --output annotated_peaks.tsv \
  --promoter-window 2000
```

This identifies:
- Promoters
- Gene bodies
- Intergenic regions
- Enhancers
- Other genomic features

## Troubleshooting

### Common Issues

1. **High mitochondrial contamination**: Check mitochondrial filtering steps
2. **Unusual fragment size distribution**: Verify Tn5 shifting and examine fragment size plots
3. **Poor peak calling**: Check for sufficient sequencing depth and proper controls

### Quality Control Metrics

Important QC metrics for ATAC-seq:
- Percentage of reads in peaks
- TSS enrichment score
- Fragment size distribution
- Library complexity

## References

- nf-core ATAC-seq pipeline: [https://nf-co.re/atacseq](https://nf-co.re/atacseq)
- ENCODE ATAC-seq guidelines: [https://www.encodeproject.org/atac-seq/](https://www.encodeproject.org/atac-seq/)
- Buenrostro et al. protocol: [https://doi.org/10.1038/nprot.2013.118](https://doi.org/10.1038/nprot.2013.118)
