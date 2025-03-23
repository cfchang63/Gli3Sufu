# Gli3-Sufu Gene Regulation Analysis 

A comprehensive toolkit for analyzing gene regulation data, including ATAC-seq, ChIP-seq, and RNA-seq.

## Overview

This repository contains scripts and workflows for analyzing various types of genomic data related to gene regulation:

- **ATAC-seq**: To identify regions of open chromatin
- **ChIP-seq**: To identify protein-DNA interactions using Gli3 pulldown
- **CUT&RUN**: To identify protein-DNA interactions using Gli3 pulldown


Each analysis type has its own pipeline with standardized inputs and outputs, making it easy to integrate results across different data types.

## Repository Structure

```
gene-regulation/
├── scripts/
│   ├── atac-seq/           # ATAC-seq analysis scripts
│   │   ├── pipeline.sh     # Main ATAC-seq pipeline script
│   │   └── utils/          # Utility scripts for ATAC-seq
│   ├── chip-seq/           # ChIP-seq analysis scripts
│   ├── rna-seq/            # RNA-seq analysis scripts
│   └── common/             # Shared utility scripts
├── config/                 # Configuration files
│   ├── genomes.config      # Reference genome configurations
│   └── defaults.config     # Default pipeline parameters
├── examples/               # Example data and usage examples
├── docs/                   # Documentation
│   ├── atac-seq.md         # ATAC-seq analysis documentation
│   ├── chip-seq.md         # ChIP-seq analysis documentation
│   └── rna-seq.md          # RNA-seq analysis documentation
└── environment.yml         # Conda environment specification
```

## Installation

### Prerequisites

- Bash (>= 4.0)
- Python (>= 3.8)
- R (>= 4.0)
- Conda/Miniconda (recommended for dependency management)

### Setting up the environment

```bash
# Clone the repository
git clone https://github.com/cfchang63/Gli3Sufu.git
cd Gli3Sufu

# Create and activate the conda environment
conda env create -f environment.yml
conda activate Gli3Sufu
```

## Usage

### ATAC-seq Analysis

The ATAC-seq pipeline processes raw FASTQ files through quality control, alignment, filtering, and peak calling:

```bash
# Basic usage with default parameters
cd scripts/atac-seq
./pipeline.sh --samples samples.txt --output results

# Advanced usage with custom parameters
./pipeline.sh --samples samples.txt --output results --genome mm10 --threads 16
```

The `samples.txt` file should contain one sample name per line. The pipeline expects paired-end FASTQ files named as `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz` in the current directory.

### ChIP-seq Analysis

[Documentation placeholder for ChIP-seq]

### CUT&RUN Analysis

[Documentation placeholder for CUT&RUN]

## Outputs

### ATAC-seq Pipeline

The ATAC-seq pipeline generates the following output directories:

- `fastqc/`: Quality control reports
- `trimmed/`: Adapter-trimmed FASTQ files
- `aligned/`: BAM files from alignment
- `filtered/`: Filtered BAM files (no mitochondria, duplicates, etc.)
- `peaks/`: Called peaks and associated files
- `bigwig/`: BigWig coverage tracks
- `differential/`: Differential peak analysis between conditions

## Data Integration

[Documentation placeholder for integrating different data types]

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the [MIT License](LICENSE).

## Citation

If you use this software in your research, please cite:

```
[Citation information placeholder]
```

## Contact

Ching-Fang Chang - [Ching-Fang.Chang@cchmc.org]

Project Link: [https://github.com/cfchang63/Gli3Sufu](https://github.com/cfchang63/Gli3Sufu)
#
