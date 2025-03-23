#!/usr/bin/env bash
#
# Gli3-Sufu Analysis Master Workflow
# 
# This script orchestrates the complete analysis workflow for
# integrating ChIP-seq, CUT&RUN, and ATAC-seq data for Gli3-Sufu
# gene regulatory network analysis.
#
# Author: Claude
# Date: March 23, 2025
#

set -e  # Exit on error
set -o pipefail  # Exit if any command in pipe fails

# CONFIGURATION
##############################################################################

# Default parameters (can be overridden via command line arguments)
GENOME="mm10"
OUTPUT_DIR="./results"
THREADS=16
SAMPLE_SHEET=""
WORKFLOW="cutandrun"
RUN_MOTIF=true
RUN_ENRICHMENT=true
GENE_LIST=""

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -o, --output DIR         Output directory (default: ${OUTPUT_DIR})"
    echo "  -g, --genome NAME        Genome name (default: ${GENOME})"
    echo "  -t, --threads NUM        Number of threads (default: ${THREADS})"
    echo "  -s, --sample-sheet FILE  Sample sheet for nf-core (required)"
    echo "  -w, --workflow NAME      nf-core workflow to use (chipseq or cutandrun, default: ${WORKFLOW})"
    echo "  --skip-motif             Skip motif analysis"
    echo "  --skip-enrichment        Skip enrichment analysis"
    echo "  --gene-list FILE         Gene list for enrichment analysis"
    echo "  -h, --help               Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
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
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -s|--sample-sheet)
            SAMPLE_SHEET="$2"
            shift 2
            ;;
        -w|--workflow)
            WORKFLOW="$2"
            shift 2
            ;;
        --skip-motif)
            RUN_MOTIF=false
            shift
            ;;
        --skip-enrichment)
            RUN_ENRICHMENT=false
            shift
            ;;
        --gene-list)
            GENE_LIST="$2"
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

# Check if sample sheet is provided
if [ -z "${SAMPLE_SHEET}" ]; then
    log "ERROR: Sample sheet is required. Use -s or --sample-sheet to specify."
    usage
fi

# Check if sample sheet exists
if [ ! -f "${SAMPLE_SHEET}" ]; then
    log "ERROR: Sample sheet file not found: ${SAMPLE_SHEET}"
    exit 1
fi

# Create main output directory
mkdir -p "${OUTPUT_DIR}"

# Step 1: Run peak intersection workflow
log "Step 1: Running peak intersection workflow..."
./peak-intersection-workflow.sh \
    --output "${OUTPUT_DIR}" \
    --genome "${GENOME}" \
    --run-nf-core \
    --sample-sheet "${SAMPLE_SHEET}" \
    --workflow "${WORKFLOW}"

# Step 2: Run motif analysis if enabled
if [ "${RUN_MOTIF}" = true ]; then
    log "Step 2: Running motif analysis..."
    MOTIF_DIR="${OUTPUT_DIR}/motif_results"
    ./motif-analysis.sh \
        --output "${MOTIF_DIR}" \
        --genome "${GENOME}" \
        --intersections "${OUTPUT_DIR}/intersections" \
        --threads "${THREADS}"
else
    log "Step 2: Skipping motif analysis..."
fi

# Step 3: Run enrichment analysis if enabled
if [ "${RUN_ENRICHMENT}" = true ]; then
    log "Step 3: Running enrichment analysis..."
    ENRICHMENT_DIR="${OUTPUT_DIR}/enrichment_results"
    
    # Process each annotation file
    for anno_file in "${OUTPUT_DIR}"/annotations/*_annotation.txt; do
        if [ -f "${anno_file}" ]; then
            basename=$(basename "${anno_file}" _annotation.txt)
            target_dir="${ENRICHMENT_DIR}/${basename}"
            
            log "  Analyzing ${basename}..."
            
            gene_list_param=""
            if [ -n "${GENE_LIST}" ]; then
                gene_list_param="--gene-list ${GENE_LIST}"
            fi
            
            python analyze_peak_enrichment.py \
                --input "${anno_file}" \
                --output "${target_dir}" \
                --genome "${GENOME}" \
                ${gene_list_param}
        fi
    done
else
    log "Step 3: Skipping enrichment analysis..."
fi

log "Analysis workflow completed successfully!"
log "Results can be found in: ${OUTPUT_DIR}"
exit 0
