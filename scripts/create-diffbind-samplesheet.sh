#!/usr/bin/env bash
#
# Create a DiffBind sample sheet for Split-DamID analysis
# 
# This script creates a sample sheet for DiffBind analysis of Split-DamID data
# based on the output from the Split-DamID workflow.
#
# Author: Ching-Fang Chang
# Date: March 23, 2025
#

set -e  # Exit on error

# CONFIGURATION
##############################################################################

# Default parameters (can be overridden via command line arguments)
SPLITDAMID_DIR="./splitdamid_results"
OUTPUT_FILE="DiffBindSampleSheet.csv"
TISSUE="NIH3T3"

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -i, --input DIR          Input directory with Split-DamID results (default: ${SPLITDAMID_DIR})"
    echo "  -o, --output FILE        Output sample sheet file (default: ${OUTPUT_FILE})"
    echo "  -t, --tissue NAME        Tissue/cell type (default: ${TISSUE})"
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
        -i|--input)
            SPLITDAMID_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -t|--tissue)
            TISSUE="$2"
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

# Check if input directory exists
if [ ! -d "${SPLITDAMID_DIR}" ]; then
    log "ERROR: Input directory ${SPLITDAMID_DIR} not found."
    exit 1
fi

# Find all sample directories in the peaks directory
log "Looking for sample peak directories in ${SPLITDAMID_DIR}/peaks"
SAMPLE_DIRS=$(find "${SPLITDAMID_DIR}/peaks" -maxdepth 1 -type d -not -path "${SPLITDAMID_DIR}/peaks" | sort)

if [ -z "${SAMPLE_DIRS}" ]; then
    log "ERROR: No sample directories found in ${SPLITDAMID_DIR}/peaks"
    exit 1
fi

# Create CSV header
log "Creating DiffBind sample sheet: ${OUTPUT_FILE}"
echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,Peaks,PeakCaller" > "${OUTPUT_FILE}"

# Process each sample
for SAMPLE_DIR in ${SAMPLE_DIRS}; do
    SAMPLE=$(basename "${SAMPLE_DIR}")
    
    # Extract condition and replicate from sample name
    if [[ "${SAMPLE}" =~ Gli3-Hand2-Dox-([0-9]+) ]]; then
        CONDITION="Gli3-Hand2-Dox"
        REPLICATE="${BASH_REMATCH[1]}"
        FACTOR="Gli3-Hand2"
    elif [[ "${SAMPLE}" =~ Gli3-Hand2-([0-9]+) ]]; then
        CONDITION="Gli3-Hand2"
        REPLICATE="${BASH_REMATCH[1]}"
        FACTOR="Gli3-Hand2"
    elif [[ "${SAMPLE}" =~ DAM-Dox-([0-9]+) ]]; then
        CONDITION="DAM-Dox"
        REPLICATE="${BASH_REMATCH[1]}"
        FACTOR="DAM"
    elif [[ "${SAMPLE}" =~ DAM-([0-9]+) ]]; then
        CONDITION="DAM"
        REPLICATE="${BASH_REMATCH[1]}"
        FACTOR="DAM"
    else
        log "WARNING: Unable to determine condition and replicate for ${SAMPLE}, skipping"
        continue
    fi
    
    # Check if BAM file exists
    BAM_FILE="${SPLITDAMID_DIR}/bam/${SAMPLE}_filtered_sorted.bam"
    if [ ! -f "${BAM_FILE}" ]; then
        log "WARNING: BAM file not found for ${SAMPLE}: ${BAM_FILE}"
        continue
    fi
    
    # Check if peak file exists
    PEAK_FILE="${SPLITDAMID_DIR}/peaks/${SAMPLE}/${SAMPLE}_broad_peaks.broadPeak"
    if [ ! -f "${PEAK_FILE}" ]; then
        # Try alternative peak file formats
        PEAK_FILE="${SPLITDAMID_DIR}/peaks/${SAMPLE}/${SAMPLE}_peaks.broadPeak"
        if [ ! -f "${PEAK_FILE}" ]; then
            PEAK_FILE="${SPLITDAMID_DIR}/peaks/${SAMPLE}/${SAMPLE}_broad_peaks.bed"
            if [ ! -f "${PEAK_FILE}" ]; then
                log "WARNING: Peak file not found for ${SAMPLE}"
                continue
            fi
        fi
    fi
    
    # Add sample to sample sheet
    echo "${SAMPLE},${TISSUE},${FACTOR},${CONDITION},${REPLICATE},${BAM_FILE},${PEAK_FILE},bed" >> "${OUTPUT_FILE}"
    log "Added sample: ${SAMPLE} (${CONDITION}, replicate ${REPLICATE})"
done

# Count samples
SAMPLE_COUNT=$(grep -v "SampleID" "${OUTPUT_FILE}" | wc -l)

if [ ${SAMPLE_COUNT} -eq 0 ]; then
    log "ERROR: No valid samples found. Check the input directory and file paths."
    exit 1
fi

log "Successfully created DiffBind sample sheet with ${SAMPLE_COUNT} samples."
log "Output file: ${OUTPUT_FILE}"

# Provide next steps
log "Next steps:"
log "  1. Verify the sample sheet to ensure all samples are correctly categorized."
log "  2. Run the differential binding analysis with:"
log "     Rscript splitdamid_diffbind_analysis.R ${OUTPUT_FILE} diffbind_results"

exit 0