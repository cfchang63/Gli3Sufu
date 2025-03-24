#!/usr/bin/env bash
#
# Peak Intersection and Annotation Workflow
# 
# This script automates the process of intersecting peak files from 
# different ChIP-seq, CUT&RUN, and ATAC-seq experiments and annotating
# the results using HOMER's annotatePeaks.pl.
#
# Author: Ching-Fang Chang
# Date: March 23, 2025
#

set -e  # Exit on error
set -o pipefail  # Exit if any command in pipe fails

# CONFIGURATION
##############################################################################

# Default parameters (can be overridden via command line arguments)
GENOME="mm10"
HOMER_PATH="annotatePeaks.pl"
GTF_FILE="${HOME}/Databases/${GENOME}/${GENOME}.ncbiRefSeq.gtf"
OUTPUT_DIR="./results"
FASTA_FILE="${HOME}/Databases/${GENOME}/${GENOME}.fa"
BLACKLIST="${HOME}/Databases/${GENOME}/${GENOME}-blacklist.v2.bed"
BOWTIE2_INDEX="${HOME}/Databases/${GENOME}"

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -o, --output DIR         Output directory (default: ${OUTPUT_DIR})"
    echo "  -g, --genome NAME        Genome name (default: ${GENOME})"
    echo "  --gtf FILE               GTF file path (default: ${GTF_FILE})"
    echo "  --fasta FILE             Genome FASTA file (default: ${FASTA_FILE})"
    echo "  --blacklist FILE         Blacklist regions BED file (default: ${BLACKLIST})"
    echo "  --bowtie2 DIR            Bowtie2 index directory (default: ${BOWTIE2_INDEX})"
    echo "  --homer PATH             Path to HOMER's annotatePeaks.pl (default: ${HOMER_PATH})"
    echo "  --run-nf-core            Run nf-core workflow (default: skip if result files exist)"
    echo "  --sample-sheet FILE      Sample sheet for nf-core workflow (CSV format)"
    echo "  --workflow NAME          nf-core workflow to use (chipseq or cutandrun)"
    echo "  -h, --help               Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    local dependencies=(bedtools sort awk)
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
        "${OUTPUT_DIR}/bed_files"
        "${OUTPUT_DIR}/intersections"
        "${OUTPUT_DIR}/annotations"
    )
    
    log "Creating output directories..."
    for dir in "${directories[@]}"; do
        mkdir -p "${dir}"
    done
}

run_nf_core_workflow() {
    local workflow="$1"
    local sample_sheet="$2"
    
    if [ ! -f "${sample_sheet}" ]; then
        log "ERROR: Sample sheet file not found: ${sample_sheet}"
        exit 1
    fi
    
    log "Checking if conda environment 'nf-core' exists..."
    if ! conda info --envs | grep -q "nf-core"; then
        log "Creating conda environment 'nf-core'..."
        conda create --yes --name nf-core python=3.12 nf-core nextflow
    fi
    
    log "Activating conda environment 'nf-core'..."
    eval "$(conda shell.bash hook)"
    conda activate nf-core
    
    log "Running nf-core/${workflow} workflow..."
    nextflow run "nf-core/${workflow}" \
        -profile singularity \
        --input "${sample_sheet}" \
        --genome "${GENOME}" \
        --bowtie2 "${BOWTIE2_INDEX}" \
        --gtf "${GTF_FILE}" \
        --blacklist "${BLACKLIST}" \
        --fasta "${FASTA_FILE}" \
        --outdir "${OUTPUT_DIR}/nf-core-${workflow}"
    
    if [ "${workflow}" == "cutandrun" ]; then
        log "Setting cutandrun-specific parameters..."
        nextflow run "nf-core/${workflow}" \
            -profile singularity \
            --input "${sample_sheet}" \
            --peakcaller seacr \
            --genome "${GENOME}" \
            --bowtie2 "${BOWTIE2_INDEX}" \
            --gtf "${GTF_FILE}" \
            --blacklist "${BLACKLIST}" \
            --fasta "${FASTA_FILE}" \
            --seacr_stringent relaxed \
            --outdir "${OUTPUT_DIR}/nf-core-${workflow}"
    fi
    
    log "nf-core/${workflow} workflow completed."
}

prepare_bed_file() {
    local input_file="$1"
    local output_file="$2"
    
    log "Preparing BED file: ${input_file} -> ${output_file}"
    
    # Extract first 3 columns (chromosome, start, end)
    cut -f1-3 "${input_file}" > "${output_file}.tmp"
    
    # Sort BED file
    sort -k1,1 -k2,2n "${output_file}.tmp" > "${output_file}"
    
    # Clean up
    rm "${output_file}.tmp"
    
    log "BED file prepared: ${output_file}"
}

intersect_bed_files() {
    local file_a="$1"
    local file_b="$2"
    local output_file="$3"
    
    log "Intersecting BED files: ${file_a} vs ${file_b} -> ${output_file}"
    
    bedtools intersect -a "${file_a}" -b "${file_b}" -wa -wb > "${output_file}"
    
    log "Intersection complete: ${output_file}"
}

annotate_peaks() {
    local bed_file="$1"
    local output_file="$2"
    local genome="$3"
    
    log "Annotating peaks: ${bed_file} -> ${output_file}"
    
    "${HOMER_PATH}" "${bed_file}" "${genome}" > "${output_file}"
    
    log "Annotation complete: ${output_file}"
}

# MAIN SCRIPT
##############################################################################

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
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
        --homer)
            HOMER_PATH="$2"
            shift 2
            ;;
        --run-nf-core)
            RUN_NF_CORE=true
            shift
            ;;
        --sample-sheet)
            SAMPLE_SHEET="$2"
            shift 2
            ;;
        --workflow)
            WORKFLOW="$2"
            shift 2
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

# Check dependencies
check_dependencies

# Create output directories
create_directories

# Run nf-core workflow if requested
if [ "${RUN_NF_CORE}" = true ]; then
    if [ -z "${WORKFLOW}" ]; then
        log "ERROR: Workflow not specified with --workflow."
        exit 1
    fi
    
    if [ -z "${SAMPLE_SHEET}" ]; then
        log "ERROR: Sample sheet not specified with --sample-sheet."
        exit 1
    fi
    
    run_nf_core_workflow "${WORKFLOW}" "${SAMPLE_SHEET}"
fi

# Prepare BED files section
log "Preparing BED files..."

# Define input files (modify paths as needed)
GLI3_INPUT="${OUTPUT_DIR}/bed_files/gli3_input.bed"
FOXF2_INPUT="${OUTPUT_DIR}/bed_files/foxf2_input.bed"
ATACSEQ_INPUT="${OUTPUT_DIR}/bed_files/atacseq_input.bed"
SATB2_INPUT="${OUTPUT_DIR}/bed_files/satb2_input.bed"

# Ask user for input files if they don't exist
if [ ! -f "${GLI3_INPUT}" ]; then
    read -p "Enter path to Gli3 peak BED file: " GLI3_PATH
    prepare_bed_file "${GLI3_PATH}" "${GLI3_INPUT}"
fi

if [ ! -f "${FOXF2_INPUT}" ]; then
    read -p "Enter path to Foxf2 peak BED file: " FOXF2_PATH
    prepare_bed_file "${FOXF2_PATH}" "${FOXF2_INPUT}"
fi

if [ ! -f "${ATACSEQ_INPUT}" ]; then
    read -p "Enter path to ATAC-seq peak BED file: " ATACSEQ_PATH
    prepare_bed_file "${ATACSEQ_PATH}" "${ATACSEQ_INPUT}"
fi

if [ ! -f "${SATB2_INPUT}" ]; then
    read -p "Enter path to SATB2 peak BED file: " SATB2_PATH
    prepare_bed_file "${SATB2_PATH}" "${SATB2_INPUT}"
fi

# Perform intersections
log "Performing peak intersections..."

# Gli3 and Foxf2
INTERSECTION_GLI3_FOXF2="${OUTPUT_DIR}/intersections/Gli3_Foxf2.bed"
intersect_bed_files "${GLI3_INPUT}" "${FOXF2_INPUT}" "${INTERSECTION_GLI3_FOXF2}"
annotate_peaks "${INTERSECTION_GLI3_FOXF2}" "${OUTPUT_DIR}/annotations/Gli3_Foxf2_annotation.txt" "${GENOME}"

# Gli3 and SATB2
INTERSECTION_GLI3_SATB2="${OUTPUT_DIR}/intersections/Gli3_Satb2.bed"
intersect_bed_files "${GLI3_INPUT}" "${SATB2_INPUT}" "${INTERSECTION_GLI3_SATB2}"
annotate_peaks "${INTERSECTION_GLI3_SATB2}" "${OUTPUT_DIR}/annotations/Gli3_Satb2_annotation.txt" "${GENOME}"

# Foxf2 and SATB2
INTERSECTION_FOXF2_SATB2="${OUTPUT_DIR}/intersections/Foxf2_Satb2.bed"
intersect_bed_files "${FOXF2_INPUT}" "${SATB2_INPUT}" "${INTERSECTION_FOXF2_SATB2}"
annotate_peaks "${INTERSECTION_FOXF2_SATB2}" "${OUTPUT_DIR}/annotations/Foxf2_Satb2_annotation.txt" "${GENOME}"

# Gli3, Foxf2, and ATAC-seq
INTERSECTION_GLI3_FOXF2_ATACSEQ="${OUTPUT_DIR}/intersections/Gli3_Foxf2_ATACseq.bed"
intersect_bed_files "${INTERSECTION_GLI3_FOXF2}" "${ATACSEQ_INPUT}" "${INTERSECTION_GLI3_FOXF2_ATACSEQ}"
annotate_peaks "${INTERSECTION_GLI3_FOXF2_ATACSEQ}" "${OUTPUT_DIR}/annotations/Gli3_Foxf2_ATACseq_annotation.txt" "${GENOME}"

# Gli3, Foxf2, ATAC-seq, and SATB2
INTERSECTION_GLI3_FOXF2_ATACSEQ_SATB2="${OUTPUT_DIR}/intersections/Gli3_Foxf2_ATACseq_Satb2.bed"
intersect_bed_files "${INTERSECTION_GLI3_FOXF2_ATACSEQ}" "${SATB2_INPUT}" "${INTERSECTION_GLI3_FOXF2_ATACSEQ_SATB2}"
annotate_peaks "${INTERSECTION_GLI3_FOXF2_ATACSEQ_SATB2}" "${OUTPUT_DIR}/annotations/Gli3_Foxf2_ATACseq_Satb2_annotation.txt" "${GENOME}"

log "Processing completed successfully."
log "Results can be found in: ${OUTPUT_DIR}/annotations/"

# Print summary of results
echo ""
echo "Summary of results:"
echo "==================="
echo "Gli3 and Foxf2 overlapping peaks: $(wc -l < ${INTERSECTION_GLI3_FOXF2})"
echo "Gli3 and SATB2 overlapping peaks: $(wc -l < ${INTERSECTION_GLI3_SATB2})"
echo "Foxf2 and SATB2 overlapping peaks: $(wc -l < ${INTERSECTION_FOXF2_SATB2})"
echo "Gli3, Foxf2, and ATAC-seq overlapping peaks: $(wc -l < ${INTERSECTION_GLI3_FOXF2_ATACSEQ})"
echo "Gli3, Foxf2, ATAC-seq, and SATB2 overlapping peaks: $(wc -l < ${INTERSECTION_GLI3_FOXF2_ATACSEQ_SATB2})"
echo ""

exit 0
