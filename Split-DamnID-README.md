# Split-DamID Analysis Workflow

A comprehensive toolkit for analyzing Split-DamID data to identify protein-DNA interactions in vivo.

## Overview

Split-DamID is a technique for mapping protein-DNA interactions in vivo without antibodies or fixation. This workflow automates the analysis of Split-DamID data, from raw sequencing reads to identifying binding sites and calculating log2 fold changes between experimental and control samples.

The workflow consists of:

1. **Read Preprocessing**: Quality control and adapter trimming
2. **Alignment**: Using Bowtie2 to align reads to the genome
3. **GATC Extension**: Extending reads to the nearest GATC sites in both directions
4. **Filtering**: Removing reads from blacklisted regions and mitochondrial DNA
5. **Coverage Analysis**: Generating coverage tracks and binned analysis
6. **Peak Calling**: Identifying enriched regions using MACS3
7. **Log2FC Calculation**: Computing log2 fold change between experimental and control samples

## Scripts

### 1. Main Workflow (`run-splitdamid-workflow.sh`)

This script orchestrates the entire Split-DamID analysis workflow:

```bash
./run-splitdamid-workflow.sh --samples samples.txt --output results_dir [OPTIONS]
```

### 2. GATC Extension Script (`extend_reads_to_GATC_revised.py`)

This Python script extends reads to the nearest GATC sites in both directions, which is essential for Split-DamID analysis:

```bash
python extend_reads_to_GATC_revised.py input.bedpe gatc_sites.bed output.bedpe
```

### 3. GATC Sites Creation (`create_gatc_sites.py`)

This script creates a BED file of all GATC restriction sites in a genome:

```bash
python create_gatc_sites.py genome.fa output.bed
```

### 4. Blacklist + mtDNA Creation (`create_blacklist_mtdna.sh`)

This script creates a combined blacklist and mitochondrial DNA exclusion file:

```bash
./create_blacklist_mtdna.sh --genome mm10 --output blacklist_mtDNA.bed
```

## Prerequisites

Before running the workflow, you need to prepare:

1. **Genome Files**:
   - Reference genome FASTA file
   - Bowtie2 index
   - GATC sites BED file (can be created using `create_gatc_sites.py`)
   - Blacklist + mtDNA exclusion file (can be created using `create_blacklist_mtdna.sh`)

2. **Sample Naming**:
   The workflow expects a specific naming pattern for control and experimental samples:
   - DAM-only controls (e.g., `DAM-1`, `DAM-2`, `DAM-3`)
   - DAM-fusion experimental samples (e.g., `Gli3-Hand2-1`, `Gli3-Hand2-2`, `Gli3-Hand2-3`)
   - Doxycycline-treated samples should include "Dox" in the name (e.g., `DAM-Dox-1`, `Gli3-Hand2-Dox-1`)

## Workflow Setup

### Step 1: Install Dependencies

```bash
# Create and activate conda environment
conda create -n splitdamid python=3.9 bowtie2 samtools bedtools deeptools r-base bioconductor-rtracklayer

# Install MACS3 in a separate environment
conda create -n MACS3 python=3.9 macs3
```

### Step 2: Create Required Files

```bash
# Create GATC sites file
python create_gatc_sites.py /path/to/genome.fa /path/to/gatc_sites.bed

# Create blacklist + mtDNA file
./create_blacklist_mtdna.sh --genome mm10 --output blacklist_mtDNA.bed
```

### Step 3: Prepare Sample File

Create a text file with one sample name per line:

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

### Step 4: Run the Workflow

```bash
./run-splitdamid-workflow.sh \
  --samples samples.txt \
  --output splitdamid_results \
  --genome mm10 \
  --fasta /path/to/mm10.fa \
  --bowtie2 /path/to/mm10_bt2index \
  --gatc /path/to/gatc_sites.bed \
  --blacklist /path/to/blacklist_mtDNA.bed \
  --extension-script /path/to/extend_reads_to_GATC_revised.py
```

## Output Structure

```
splitdamid_results/
├── bam/                   # BAM files (aligned, extended, filtered)
├── bed/                   # BED files (BEDPE format)
├── bigwig/                # BigWig coverage tracks
├── bedgraph/              # BedGraph coverage files
├── peaks/                 # Peak calls from MACS3
│   ├── Sample1/
│   └── Sample2/
├── analysis/              # Binned coverage analysis
├── log2fc/                # Log2 fold change calculations and tracks
└── genome_bins_75bp.bed   # Genome bins for coverage analysis
```

## Integrating with Other Data Types

The results from this Split-DamID analysis can be integrated with ChIP-seq, CUT&RUN, and ATAC-seq data using the peak-intersection-workflow.sh script included in the package.

To integrate Split-DamID peaks with other data types:

```bash
./peak-intersection-workflow.sh \
  --output integrated_results \
  --genome mm10 \
  --splitdamid-peaks /path/to/splitdamid_results/peaks \
  --chipseq-peaks /path/to/chipseq_results/peaks \
  --atacseq-peaks /path/to/atacseq_results/peaks
```

## Advanced Usage

### Customizing Sample Pattern Matching

You can customize how the workflow identifies control and experimental samples:

```bash
./run-splitdamid-workflow.sh \
  --control-pattern "DAM-Control" \
  --treat-pattern "DAM-Gli3" \
  --dox-pattern "DOX"
```

### Changing Bin Size

Change the bin size for coverage analysis:

```bash
./run-splitdamid-workflow.sh --bin-size 100
```

## Troubleshooting

### GATC Extension Issues

If the GATC extension step fails, check:
- The GATC sites BED file format is correct
- The extension script is using the correct Python libraries
- The input BEDPE file has the correct format

### Log2FC Calculation Issues

For log2FC calculation issues:
- Ensure R with the required packages (GenomicRanges, rtracklayer) is installed
- Check that control and experimental samples are correctly identified by pattern matching

## References

1. Kind J, et al. (2015). Genome-wide maps of nuclear lamina interactions in single human cells. Cell, 163(1), 134-147.
2. Pindyurin AV, et al. (2016). The large fraction of heterochromatin in Drosophila neurons is bound by both B-type lamin and HP1a. Epigenetics & Chromatin, 9, 22.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
