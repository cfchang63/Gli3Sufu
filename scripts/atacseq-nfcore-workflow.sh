#!/usr/bin/env bash
#
# ATAC-seq Analysis Workflow
# 
# This script automates the process of running nf-core/atacseq pipeline
# and integrating the results with ChIP-seq and CUT&RUN data.
#
# Author: Ching-Fang Chang
# Date: March 23, 2025
#

# Usage:
#  ./run-atacseq-workflow.sh \
#  --sample-sheet atacseq_samples.csv \
#  --output atacseq_results \
#  --genome mm10 \
#  --fasta /path/to/Databases/mm10/mm10.fa \
#  --gtf /path/to/Databases/mm10/mm10.ncbiRefSeq.gtf \
#  --blacklist /path/to/Databases/mm10/mm10-blacklist.v2.bed \
#  --chip-peaks /path/to/chipseq/bed_files \
#  --cutrun-peaks /path/to/cutrun/bed_files



set -e  # Exit on error
set -o pipefail  # Exit if any command in pipe fails

# CONFIGURATION
##############################################################################

# Default parameters (can be overridden via command line arguments)
GENOME="mm10"
BOWTIE2_INDEX="${HOME}/Databases/${GENOME}"
GTF_FILE="${HOME}/Databases/${GENOME}/${GENOME}.ncbiRefSeq.gtf"
FASTA_FILE="${HOME}/Databases/${GENOME}/${GENOME}.fa"
BLACKLIST="${HOME}/Databases/${GENOME}/${GENOME}-blacklist.v2.bed"
OUTPUT_DIR="./atacseq_results"
SAMPLE_SHEET=""
THREADS=16
HOMER_PATH="annotatePeaks.pl"
RUN_INTEGRATION=true

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -s, --sample-sheet FILE  Sample sheet for nf-core/atacseq (required)"
    echo "  -o, --output DIR         Output directory (default: ${OUTPUT_DIR})"
    echo "  -g, --genome NAME        Genome name (default: ${GENOME})"
    echo "  --gtf FILE               GTF file path (default: ${GTF_FILE})"
    echo "  --fasta FILE             Genome FASTA file (default: ${FASTA_FILE})"
    echo "  --blacklist FILE         Blacklist regions BED file (default: ${BLACKLIST})"
    echo "  --bowtie2 DIR            Bowtie2 index directory (default: ${BOWTIE2_INDEX})"
    echo "  -t, --threads NUM        Number of threads (default: ${THREADS})"
    echo "  --homer PATH             Path to HOMER's annotatePeaks.pl (default: ${HOMER_PATH})"
    echo "  --chip-peaks DIR         Directory with ChIP-seq peak files for integration"
    echo "  --cutrun-peaks DIR       Directory with CUT&RUN peak files for integration"
    echo "  --skip-integration       Skip integration with ChIP-seq and CUT&RUN data"
    echo "  -h, --help               Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    local dependencies=(nextflow bedtools sort)
    local missing=()
    
    log "Checking dependencies..."
    for cmd in "${dependencies[@]}"; do
        if ! command -v "${cmd}" &> /dev/null; then
            missing+=("${cmd}")
        fi
    done
    
    # Check if HOMER's annotatePeaks.pl is available
    if ! command -v "${HOMER_PATH}" &> /dev/null; then
        log "WARNING: HOMER's annotatePeaks.pl not found in PATH. Will use specified path: ${HOMER_PATH}"
        if [ ! -x "${HOMER_PATH}" ]; then
            missing+=("annotatePeaks.pl")
        fi
    fi
    
    if [ ${#missing[@]} -gt 0 ]; then
        log "ERROR: Missing dependencies: ${missing[*]}"
        log "Please install missing dependencies and try again."
        exit 1
    fi
    
    log "All dependencies found."
}

create_directories() {
    local directories=(
        "${OUTPUT_DIR}"
        "${OUTPUT_DIR}/integrations"
        "${OUTPUT_DIR}/annotations"
    )
    
    log "Creating output directories..."
    for dir in "${directories[@]}"; do
        mkdir -p "${dir}"
    done
    
    log "Output directories created."
}

check_conda_env() {
    log "Checking if conda environment 'nf-core' exists..."
    if ! conda info --envs | grep -q "nf-core"; then
        log "Creating conda environment 'nf-core'..."
        conda create --yes --name nf-core python=3.12 nf-core nextflow
    fi
    
    log "Activating conda environment 'nf-core'..."
    eval "$(conda shell.bash hook)"
    conda activate nf-core
}

run_atacseq_pipeline() {
    local sample_sheet="$1"
    
    log "Running nf-core/atacseq pipeline..."
    nextflow run nf-core/atacseq \
        -profile singularity \
        --input "${sample_sheet}" \
        --genome "${GENOME}" \
        --fasta "${FASTA_FILE}" \
        --gtf "${GTF_FILE}" \
        --blacklist "${BLACKLIST}" \
        --outdir "${OUTPUT_DIR}/nf-core-atacseq" \
        --max_cpus "${THREADS}" \
        --max_memory "120.GB" \
        --macs_gsize 2.7e9 \
        --skip_diff_analysis
    
    log "nf-core/atacseq pipeline completed successfully."
}

prepare_atacseq_peaks() {
    log "Preparing ATAC-seq peaks for integration..."
    
    # Locate consensus peaks from nf-core/atacseq output
    local consensus_peaks="${OUTPUT_DIR}/nf-core-atacseq/macs/consensus/consensus_peaks.mLb.clN.bed"
    
    if [ ! -f "${consensus_peaks}" ]; then
        log "ERROR: Could not find consensus peaks file at ${consensus_peaks}"
        log "Checking for alternative locations..."
        
        # Try to find the consensus peaks file
        local alternative_peaks=$(find "${OUTPUT_DIR}/nf-core-atacseq" -name "consensus_peaks*.bed" | head -n 1)
        
        if [ -z "${alternative_peaks}" ]; then
            log "ERROR: Could not find any consensus peaks file. Please check the pipeline output."
            exit 1
        else
            consensus_peaks="${alternative_peaks}"
            log "Found consensus peaks at ${consensus_peaks}"
        fi
    fi
    
    # Extract coordinates
    cut -f1-3 "${consensus_peaks}" > "${OUTPUT_DIR}/atacseq_peaks.bed"
    
    # Sort BED file
    sort -k1,1 -k2,2n "${OUTPUT_DIR}/atacseq_peaks.bed" > "${OUTPUT_DIR}/atacseq_peaks_sorted.bed"
    
    log "ATAC-seq peaks prepared: ${OUTPUT_DIR}/atacseq_peaks_sorted.bed"
}

integrate_with_chipseq_cutrun() {
    local atacseq_peaks="${OUTPUT_DIR}/atacseq_peaks_sorted.bed"
    
    # Process ChIP-seq peaks if provided
    if [ -n "${CHIP_PEAKS_DIR}" ]; then
        log "Integrating ATAC-seq with ChIP-seq peaks..."
        
        # Find all BED files in the ChIP-seq directory
        chip_files=("${CHIP_PEAKS_DIR}"/*.bed)
        
        if [ ${#chip_files[@]} -eq 0 ]; then
            log "WARNING: No BED files found in ${CHIP_PEAKS_DIR}"
        else
            log "Found ${#chip_files[@]} ChIP-seq peak files"
            
            # Process each ChIP-seq file
            for chip_file in "${chip_files[@]}"; do
                local name=$(basename "${chip_file}" .bed)
                
                log "Processing ${name}..."
                
                # Prepare sorted BED file
                sort -k1,1 -k2,2n "${chip_file}" > "${OUTPUT_DIR}/integrations/${name}_sorted.bed"
                
                # Intersect with ATAC-seq peaks
                bedtools intersect -a "${OUTPUT_DIR}/integrations/${name}_sorted.bed" -b "${atacseq_peaks}" -wa -wb \
                    > "${OUTPUT_DIR}/integrations/${name}_atacseq_overlap.bed"
                
                # Annotate intersections
                "${HOMER_PATH}" "${OUTPUT_DIR}/integrations/${name}_atacseq_overlap.bed" "${GENOME}" \
                    > "${OUTPUT_DIR}/annotations/${name}_atacseq_annotation.txt"
                
                log "Integration completed for ${name}"
            done
        fi
    fi
    
    # Process CUT&RUN peaks if provided
    if [ -n "${CUTRUN_PEAKS_DIR}" ]; then
        log "Integrating ATAC-seq with CUT&RUN peaks..."
        
        # Find all BED files in the CUT&RUN directory
        cutrun_files=("${CUTRUN_PEAKS_DIR}"/*.bed)
        
        if [ ${#cutrun_files[@]} -eq 0 ]; then
            log "WARNING: No BED files found in ${CUTRUN_PEAKS_DIR}"
        else
            log "Found ${#cutrun_files[@]} CUT&RUN peak files"
            
            # Process each CUT&RUN file
            for cutrun_file in "${cutrun_files[@]}"; do
                local name=$(basename "${cutrun_file}" .bed)
                
                log "Processing ${name}..."
                
                # Prepare sorted BED file
                sort -k1,1 -k2,2n "${cutrun_file}" > "${OUTPUT_DIR}/integrations/${name}_sorted.bed"
                
                # Intersect with ATAC-seq peaks
                bedtools intersect -a "${OUTPUT_DIR}/integrations/${name}_sorted.bed" -b "${atacseq_peaks}" -wa -wb \
                    > "${OUTPUT_DIR}/integrations/${name}_atacseq_overlap.bed"
                
                # Annotate intersections
                "${HOMER_PATH}" "${OUTPUT_DIR}/integrations/${name}_atacseq_overlap.bed" "${GENOME}" \
                    > "${OUTPUT_DIR}/annotations/${name}_atacseq_annotation.txt"
                
                log "Integration completed for ${name}"
            done
        fi
    fi
    
    # Annotate ATAC-seq peaks themselves
    log "Annotating ATAC-seq peaks..."
    "${HOMER_PATH}" "${atacseq_peaks}" "${GENOME}" > "${OUTPUT_DIR}/annotations/atacseq_annotation.txt"
    
    log "Integration and annotation completed."
}

create_sample_sheet_template() {
    local template="${OUTPUT_DIR}/atacseq_sample_template.csv"
    
    log "Creating sample sheet template at ${template}..."
    
    cat > "${template}" << EOF
sample,fastq_1,fastq_2,replicate,single_end
SAMPLE1,/path/to/SAMPLE1_R1.fastq.gz,/path/to/SAMPLE1_R2.fastq.gz,1,0
SAMPLE2,/path/to/SAMPLE2_R1.fastq.gz,/path/to/SAMPLE2_R2.fastq.gz,1,0
SAMPLE3,/path/to/SAMPLE3_R1.fastq.gz,/path/to/SAMPLE3_R2.fastq.gz,2,0
SAMPLE4,/path/to/SAMPLE4_R1.fastq.gz,/path/to/SAMPLE4_R2.fastq.gz,2,0
EOF
    
    log "Sample sheet template created. Edit this file with your sample information."
}

# MAIN SCRIPT
##############################################################################

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--sample-sheet)
            SAMPLE_SHEET="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        --gtf)
            GTF_FILE="$2"
            shift 2
            ;;
        --fasta)
            FASTA_FILE="$2"
            shift 2
            ;;
        --blacklist)
            BLACKLIST="$2"
            shift 2
            ;;
        --bowtie2)
            BOWTIE2_INDEX="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --homer)
            HOMER_PATH="$2"
            shift 2
            ;;
        --chip-peaks)
            CHIP_PEAKS_DIR="$2"
            shift 2
            ;;
        --cutrun-peaks)
            CUTRUN_PEAKS_DIR="$2"
            shift 2
            ;;
        --skip-integration)
            RUN_INTEGRATION=false
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if sample sheet is provided
if [ -z "${SAMPLE_SHEET}" ]; then
    log "ERROR: Sample sheet is required. Use -s or --sample-sheet to specify."
    log "Creating a sample sheet template for your reference..."
    create_sample_sheet_template
    log "Please fill in the template and run the script again with -s option."
    exit 1
fi

# Check if sample sheet exists
if [ ! -f "${SAMPLE_SHEET}" ]; then
    log "ERROR: Sample sheet file not found: ${SAMPLE_SHEET}"
    exit 1
fi

# Check dependencies
check_dependencies

# Create output directories
create_directories

# Check and activate conda environment
check_conda_env

# Run nf-core/atacseq pipeline
run_atacseq_pipeline "${SAMPLE_SHEET}"

# Prepare ATAC-seq peaks for integration
prepare_atacseq_peaks

# Integrate with ChIP-seq and CUT&RUN data if requested
if [ "${RUN_INTEGRATION}" = true ]; then
    if [ -z "${CHIP_PEAKS_DIR}" ] && [ -z "${CUTRUN_PEAKS_DIR}" ]; then
        log "WARNING: Neither ChIP-seq nor CUT&RUN peak directories specified."
        log "Skipping integration. Use --chip-peaks and/or --cutrun-peaks to specify."
    else
        integrate_with_chipseq_cutrun
    fi
fi

log "ATAC-seq analysis workflow completed successfully!"
log "Results can be found in: ${OUTPUT_DIR}"

# Print summary
echo ""
echo "Summary of results:"
echo "==================="
echo "ATAC-seq pipeline output: ${OUTPUT_DIR}/nf-core-atacseq/"
echo "ATAC-seq peaks: ${OUTPUT_DIR}/atacseq_peaks_sorted.bed"
echo "ATAC-seq peak annotations: ${OUTPUT_DIR}/annotations/atacseq_annotation.txt"

if [ "${RUN_INTEGRATION}" = true ]; then
    echo "Integration results: ${OUTPUT_DIR}/integrations/"
    echo "Annotated integrations: ${OUTPUT_DIR}/annotations/"
fi

echo ""
log "Done."
exit 0
