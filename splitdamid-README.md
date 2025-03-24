# Split-DamID Analysis Workflow

This document provides detailed information about the Split-DamID analysis workflow implemented in this repository.

## Overview

Split-DamID is a technique for mapping protein-DNA interactions in vivo without antibodies or fixation. This method uses a split version of the Dam methyltransferase that becomes active only when reconstituted by the interaction of two proteins of interest. The workflow processes raw Split-DamID data through alignment, extension to GATC sites, and differential analysis to identify genome-wide DNA binding sites.

## Workflow Steps

The Split-DamID workflow consists of the following steps:

1. **Read Preprocessing**: Quality control and adapter trimming
2. **Alignment**: Using Bowtie2 to align reads to the genome
3. **GATC Extension**: Extending reads to the nearest GATC sites in both directions
4. **Filtering**: Removing reads from blacklisted regions and mitochondrial DNA
5. **Coverage Analysis**: Generating coverage tracks and binned analysis
6. **Peak Calling**: Identifying enriched regions using MACS3
7. **Log2FC Calculation**: Computing log2 fold change between experimental and control samples
8. **Differential Analysis**: Statistical analysis of differential binding

## Required Input

### Sample Naming Convention

The workflow expects a specific naming pattern for control and experimental samples:
- DAM-only controls (e.g., `DAM-1`, `DAM-2`, `DAM-3`)
- DAM-fusion experimental samples (e.g., `Gli3-Hand2-1`, `Gli3-Hand2-2`, `Gli3-Hand2-3`)
- Doxycycline-treated samples should include "Dox" in the name (e.g., `DAM-Dox-1`, `Gli3-Hand2-Dox-1`)

### Sample File Format

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

### Required Files

Before running Split-DamID analysis, you need to prepare:

1. GATC sites file (containing all GATC restriction sites in the genome)
2. Blacklist + mtDNA exclusion file

## Running the Workflow

### Preparing Required Files

1. Create GATC sites file:
```bash
python scripts/create-gatc-sites.py /path/to/genome.fa /path/to/gatc_sites.bed
```

2. Create blacklist + mtDNA file:
```bash
./scripts/create-blacklist-mtdna.sh --genome mm10 --output blacklist_mtDNA.bed
```

### Basic Usage

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

### Full Options

```bash
./scripts/splitdamid-workflow.sh \
  --samples samples/splitdamid_samples.txt \
  --output splitdamid_results \
  --genome mm10 \
  --fasta /path/to/mm10.fa \
  --bowtie2 /path/to/mm10_bt2index \
  --gatc /path/to/gatc_sites.bed \
  --blacklist /path/to/blacklist_mtDNA.bed \
  --extension-script scripts/gatc-extension-script.py \
  --bin-size 75 \
  --threads 24 \
  --genome-size mm \
  --control-pattern "DAM-[0-9]" \
  --treat-pattern "Gli3-Hand2-[0-9]" \
  --dox-pattern "Dox"
```

## Output Files

The workflow generates the following output directories:

- `splitdamid_results/`: Main output directory
  - `bam/`: BAM files (aligned, extended, filtered)
  - `bed/`: BED files (BEDPE format)
  - `bigwig/`: BigWig coverage tracks
  - `bedgraph/`: BedGraph coverage files
  - `peaks/`: Peak calls from MACS3
  - `analysis/`: Binned coverage analysis
  - `log2fc/`: Log2 fold change calculations and tracks
  - `merged_peaks/`: Merged peaks across samples
  - `counts/`: Read count matrices

### Key Output Files

- `*_filtered_sorted.bam`: Filtered and sorted BAM files
- `*_extended.bed`: Reads extended to GATC sites
- `*.bw`: BigWig coverage tracks
- `*_broad_peaks.broadPeak`: Peak calls from MACS3
- `*_bin_counts.bed`: Binned genomic coverage
- `*_log2fc.txt` and `*_log2fc.bw`: Log2 fold change calculations
- `unified_peaks.bed`: Merged peaks across all samples
- `combined_counts.tsv`: Combined count matrix for differential analysis

## Split-DamID-Specific Processing

### GATC Extension

A critical step in Split-DamID analysis is extending reads to the nearest GATC sites:

```bash
python scripts/gatc-extension-script.py \
  input.bedpe gatc_sites.bed output.bedpe
```

This extends each end of a sequenced fragment to the nearest GATC site in the appropriate direction, as Dam methyltransferase only methylates GATC sites.

### Binned Analysis

The genome is divided into bins (default 75bp) for quantitative analysis:

```bash
bedtools coverage -a genome_bins_75bp.bed \
  -b sample_filtered_sorted.bam -counts > \
  sample_bin_counts.bed
```

### Log2 Fold Change Calculation

Log2 fold changes are calculated between experimental and control samples:

```bash
# In R
log2FC = log2((treat_count + 1) / (control_count + 1))
```

## Differential Binding Analysis

After processing raw data, run differential binding analysis:

```bash
./scripts/splitdamid-binning-workflow.sh \
  --input splitdamid_results \
  --output splitdamid_diffbind \
  --bin-script scripts/prepare-splitdamid-bins.py \
  --diffbind-script scripts/splitdamid-diff-analysis.r
```

This performs:
1. Preparation of binned count matrix
2. DESeq2 differential analysis
3. Annotation of differential bins
4. GO and KEGG pathway enrichment
5. Visualization (PCA, heatmaps, volcano plots)

### Complete Analysis Pipeline

For a complete analysis from raw data to statistical results:

```bash
./scripts/splitdamid-analysis-wrapper.sh \
  --raw-data raw_data \
  --output splitdamid_results \
  --final-output splitdamid_analysis \
  --genome mm10 \
  --samples samples.txt
```

## Integration with Other Data Types

Split-DamID results can be integrated with other genomic data:

```bash
./scripts/peak-intersection-workflow.sh \
  --output integration_results \
  --genome mm10 \
  --splitdamid-peaks splitdamid_results/peaks \
  --chipseq-peaks chipseq_results/peaks \
  --atacseq-peaks atacseq_results/peaks
```

## Troubleshooting

### Common Issues

1. **GATC extension failures**: Check if GATC sites file format is correct
2. **Low number of extended reads**: Verify sample preparation and Dam expression
3. **High background in DAM-only controls**: Check for leaky expression or other technical issues

### Key Metrics to Check

- Percentage of reads successfully extended to GATC sites
- Correlation between replicates
- Signal-to-noise ratio between experimental and control samples
- Number of significant differential binding regions

## References

1. Kind J, et al. (2015). Genome-wide maps of nuclear lamina interactions in single human cells. Cell, 163(1), 134-147.
2. Pindyurin AV, et al. (2016). The large fraction of heterochromatin in Drosophila neurons is bound by both B-type lamin and HP1a. Epigenetics & Chromatin, 9, 22.
3. van Schaik T, et al. (2020). Split-DamID: Measuring chromatin occupancy of interacting proteins. Genome Res. 30(7):1089-1101.
