# Gli3-Sufu Peak Analysis Workflow

A comprehensive toolkit for analyzing ChIP-seq, CUT&RUN, and ATAC-seq peaks with a focus on the Gli3-Sufu gene regulatory network.

## Overview

This repository contains scripts and workflows for integrating and analyzing multiple types of genomic data:

- **ChIP-seq**: To identify Gli3, Foxf2, and SATB2 binding sites
- **CUT&RUN**: To identify protein-DNA interactions using Gli3 antibodies
- **ATAC-seq**: To identify regions of open chromatin

The scripts in this repository automate the process of intersecting peaks from different data types, annotating the intersections, and performing enrichment analysis.

## Scripts

### 1. Peak Intersection Workflow (`peak-intersection-workflow.sh`)

This script automates the process of:
- Running nf-core pipelines for ChIP-seq or CUT&RUN analysis (optional)
- Preparing BED files for intersection
- Intersecting peaks from different experiments
- Annotating the intersected peaks using HOMER

#### Usage:

```bash
./peak-intersection-workflow.sh --output results_dir --genome mm10 [OPTIONS]
```

Options:
- `-o, --output DIR`: Output directory
- `-g, --genome NAME`: Genome name (default: mm10)
- `--gtf FILE`: GTF file path
- `--fasta FILE`: Genome FASTA file
- `--blacklist FILE`: Blacklist regions BED file
- `--bowtie2 DIR`: Bowtie2 index directory
- `--homer PATH`: Path to HOMER's annotatePeaks.pl
- `--run-nf-core`: Run nf-core workflow
- `--sample-sheet FILE`: Sample sheet for nf-core workflow
- `--workflow NAME`: nf-core workflow to use (chipseq or cutandrun)

### 2. Motif Analysis (`motif-analysis.sh`)

This script performs motif analysis on the intersection peaks using HOMER's findMotifsGenome.pl.

#### Usage:

```bash
./motif-analysis.sh --output motif_results --intersections ./results/intersections [OPTIONS]
```

Options:
- `-o, --output DIR`: Output directory
- `-g, --genome NAME`: Genome name
- `-i, --intersections DIR`: Directory with intersection BED files
- `-s, --size NUM`: Size of region for motif analysis
- `-m, --motif-len LIST`: Comma-separated list of motif lengths
- `-t, --threads NUM`: Number of threads to use

### 3. Peak Enrichment Analysis (`analyze_peak_enrichment.py`)

This Python script analyzes the enrichment of peaks in genomic features and performs statistical analysis on gene enrichment.

#### Usage:

```bash
python analyze_peak_enrichment.py --input annotated_peaks.txt --output results_dir [OPTIONS]
```

Options:
- `--input, -i`: Input annotated peak file from HOMER
- `--output, -o`: Output directory for results
- `--genome, -g`: Genome assembly (default: mm10)
- `--promoter-window`: Window size around TSS (default: 2000)
- `--background, -b`: Background peak file for enrichment analysis
- `--gene-list, -l`: List of genes of interest

## Workflow Example

Here's a complete workflow example:

```bash
# 1. Run the peak intersection workflow
./peak-intersection-workflow.sh --output ./results --genome mm10 --run-nf-core \
    --sample-sheet cutnrun_samples.csv --workflow cutandrun

# 2. Run motif analysis on the intersections
./motif-analysis.sh --output ./motif_results --intersections ./results/intersections

# 3. Analyze peak enrichment
python analyze_peak_enrichment.py --input ./results/annotations/Gli3_Foxf2_ATACseq_annotation.txt \
    --output ./enrichment_results --background ./results/annotations/Gli3_annotation.txt
```

## Sample Sheet Format

For nf-core workflows, the sample sheet should be in CSV format:

```csv
sample,fastq_1,fastq_2,antibody,control
WT_Input,/path/to/WT_Input.fastq.gz,,Input,
WT_Gli3_1,/path/to/WT_Gli3_1.fastq.gz,,Gli3,WT_Input
```

## Dependencies

- Bash (>= 4.0)
- Python (>= 3.8)
- R (>= 4.0)
- Conda (for nf-core workflows)
- Nextflow
- HOMER
- BEDTools
- Pandas, NumPy, Matplotlib, Seaborn (Python packages)

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/Gli3-Sufu-Analysis.git
cd Gli3-Sufu-Analysis

# Make scripts executable
chmod +x *.sh

# Install Python dependencies
pip install pandas numpy matplotlib seaborn scipy statsmodels pyranges genomepy

# For nf-core workflows
conda create --name nf-core python=3.12 nf-core nextflow
conda activate nf-core
```

## Output Structure

```
results/
├── bed_files/               # Prepared BED files
├── intersections/           # Intersection BED files
│   ├── Gli3_Foxf2.bed
│   ├── Gli3_Satb2.bed
│   └── ...
├── annotations/             # Annotated peaks
│   ├── Gli3_Foxf2_annotation.txt
│   ├── Gli3_Satb2_annotation.txt
│   └── ...
└── nf-core-*/               # nf-core workflow output (if run)

motif_results/               # Motif analysis results
├── Gli3_Foxf2/
│   ├── homerResults.html
│   └── ...
└── ...

enrichment_results/          # Peak enrichment analysis
├── top_genes.tsv
├── top_genes.png
├── genomic_distribution.png
└── ...
```

## References

- nf-core ChIP-seq: https://nf-co.re/chipseq
- nf-core CUT&RUN: https://nf-co.re/cutandrun
- HOMER: http://homer.ucsd.edu/homer/
- BEDTools: https://bedtools.readthedocs.io/

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributors

- Your Name - [your.email@example.com]

## Acknowledgments

- The nf-core community for their excellent pipelines
- The HOMER team for their annotation and motif analysis tools
