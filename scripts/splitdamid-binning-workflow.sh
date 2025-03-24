#!/usr/bin/env bash
#
# Split-DamID Binning and Differential Analysis Workflow
# 
# This script automates the process of preparing binned count data from
# Split-DamID workflow results for differential binding analysis.
#
# Author: Ching-Fang Chang
# Date: March 23, 2025
#

set -e  # Exit on error
set -o pipefail  # Exit if any command in pipe fails

# CONFIGURATION
##############################################################################

# Default parameters (can be overridden via command line arguments)
SPLITDAMID_DIR="./splitdamid_results"
OUTPUT_DIR="./splitdamid_diffbind"
DIFFBIND_SCRIPT="splitdamid_diffbind_analysis.R"
SAMPLE_INFO=""
BIN_SCRIPT="prepare_splitdamid_bins.py"
CUSTOM_CONTRAST=""

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -i, --input DIR          Input directory with Split-DamID results (default: ${SPLITDAMID_DIR})"
    echo "  -o, --output DIR         Output directory for differential analysis (default: ${OUTPUT_DIR})"
    echo "  --sample-info FILE       Sample information file (optional, will be created if not provided)"
    echo "  --bin-script FILE        Path to binning script (default: ${BIN_SCRIPT})"
    echo "  --diffbind-script FILE   Path to DiffBind analysis R script (default: ${DIFFBIND_SCRIPT})"
    echo "  --contrast NUM           Contrast number to analyze (optional)"
    echo "  -h, --help               Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    local dependencies=(python3 Rscript)
    local missing=()
    
    log "Checking dependencies..."
    for cmd in "${dependencies[@]}"; do
        if ! command -v "${cmd}" &> /dev/null; then
            missing+=("${cmd}")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        log "ERROR: Missing dependencies: ${missing[*]}"
        log "Please install missing dependencies and try again."
        exit 1
    fi
    
    # Check for required Python script
    if [ ! -f "${BIN_SCRIPT}" ]; then
        log "ERROR: Binning script not found at ${BIN_SCRIPT}"
        log "Please specify the correct path with --bin-script"
        exit 1
    fi
    
    # Check for required R script
    if [ ! -f "${DIFFBIND_SCRIPT}" ]; then
        log "ERROR: DiffBind analysis script not found at ${DIFFBIND_SCRIPT}"
        log "Please specify the correct path with --diffbind-script"
        exit 1
    fi
    
    log "All dependencies found."
}

create_directories() {
    mkdir -p "${OUTPUT_DIR}"
    log "Created output directory: ${OUTPUT_DIR}"
}

prepare_binned_counts() {
    log "Preparing binned count matrix..."
    
    local count_matrix="${OUTPUT_DIR}/binned_count_matrix.csv"
    local sample_info_arg=""
    
    if [ -n "${SAMPLE_INFO}" ]; then
        sample_info_arg="--sample-info ${SAMPLE_INFO}"
    else
        SAMPLE_INFO="${OUTPUT_DIR}/sample_info.csv"
        sample_info_arg="--sample-info ${SAMPLE_INFO}"
    fi
    
    # Run binning script
    python3 "${BIN_SCRIPT}" \
        --input "${SPLITDAMID_DIR}" \
        --output "${count_matrix}" \
        ${sample_info_arg} \
        --bin-dir analysis \
        --bin-pattern "*_bin_counts.bed"
    
    if [ ! -f "${count_matrix}" ]; then
        log "ERROR: Failed to create count matrix"
        exit 1
    fi
    
    if [ ! -f "${SAMPLE_INFO}" ]; then
        log "ERROR: Failed to create sample info file"
        exit 1
    fi
    
    log "Binned count matrix created: ${count_matrix}"
    log "Sample info file: ${SAMPLE_INFO}"
    
    echo "${count_matrix}"  # Return path to count matrix
}

run_diffbind_analysis() {
    local count_matrix="$1"
    local contrast_arg=""
    
    if [ -n "${CUSTOM_CONTRAST}" ]; then
        contrast_arg="${CUSTOM_CONTRAST}"
    fi
    
    log "Running differential binding analysis..."
    
    # Run R script
    Rscript "${DIFFBIND_SCRIPT}" \
        "${count_matrix}" \
        "${OUTPUT_DIR}" \
        ${contrast_arg}
    
    log "Differential binding analysis complete. Results in ${OUTPUT_DIR}"
}

verify_bin_counts() {
    log "Verifying bin count files in the Split-DamID results..."
    
    # Check if bin counts exist
    bin_files=$(find "${SPLITDAMID_DIR}/analysis" -name "*_bin_counts.bed" | wc -l)
    
    if [ "${bin_files}" -eq 0 ]; then
        log "ERROR: No bin count files found in ${SPLITDAMID_DIR}/analysis"
        log "Please run the Split-DamID workflow with the binning step first."
        exit 1
    fi
    
    log "Found ${bin_files} bin count files for processing."
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
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --sample-info)
            SAMPLE_INFO="$2"
            shift 2
            ;;
        --bin-script)
            BIN_SCRIPT="$2"
            shift 2
            ;;
        --diffbind-script)
            DIFFBIND_SCRIPT="$2"
            shift 2
            ;;
        --contrast)
            CUSTOM_CONTRAST="$2"
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

# Check if Split-DamID directory exists
if [ ! -d "${SPLITDAMID_DIR}" ]; then
    log "ERROR: Split-DamID results directory not found: ${SPLITDAMID_DIR}"
    exit 1
fi

# Check dependencies
check_dependencies

# Create output directory
create_directories

# Verify bin count files exist
verify_bin_counts

# Prepare binned counts
count_matrix=$(prepare_binned_counts)

# Run differential binding analysis
run_diffbind_analysis "${count_matrix}"

log "Split-DamID binning and differential analysis workflow completed successfully!"
log "Results can be found in: ${OUTPUT_DIR}"

exit 0
