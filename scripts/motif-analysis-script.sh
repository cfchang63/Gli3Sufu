#!/usr/bin/env bash
#
# Peak Motif Analysis Workflow
# 
# This script performs motif analysis on the intersection peaks
# using HOMER's findMotifsGenome.pl
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
HOMER_PATH="findMotifsGenome.pl"
OUTPUT_DIR="./motif_results"
INTERSECTION_DIR="./results/intersections"
SIZE=200
MOTIF_LENGTH="8,10,12"
THREADS=16

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -o, --output DIR         Output directory (default: ${OUTPUT_DIR})"
    echo "  -g, --genome NAME        Genome name (default: ${GENOME})"
    echo "  -i, --intersections DIR  Directory with intersection BED files (default: ${INTERSECTION_DIR})"
    echo "  -s, --size NUM           Size of region for motif analysis (default: ${SIZE})"
    echo "  -m, --motif-len LIST     Comma-separated list of motif lengths (default: ${MOTIF_LENGTH})"
    echo "  -t, --threads NUM        Number of threads to use (default: ${THREADS})"
    echo "  --homer PATH             Path to HOMER's findMotifsGenome.pl (default: ${HOMER_PATH})"
    echo "  -h, --help               Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    # Check if HOMER's findMotifsGenome.pl is available
    if ! command -v "${HOMER_PATH}" &> /dev/null; then
        log "WARNING: HOMER's findMotifsGenome.pl not found in PATH. Will use specified path: ${HOMER_PATH}"
        if [ ! -x "${HOMER_PATH}" ]; then
            log "ERROR: HOMER's findMotifsGenome.pl not found at specified path: ${HOMER_PATH}"
            exit 1
        fi
    fi
    
    log "All dependencies found."
}

create_directories() {
    mkdir -p "${OUTPUT_DIR}"
    log "Created output directory: ${OUTPUT_DIR}"
}

run_motif_analysis() {
    local bed_file="$1"
    local output_dir="$2"
    local name=$(basename "${bed_file}" .bed)
    
    log "Running motif analysis on ${bed_file}..."
    
    mkdir -p "${output_dir}/${name}"
    
    "${HOMER_PATH}" "${bed_file}" "${GENOME}" "${output_dir}/${name}" \
        -size "${SIZE}" \
        -len "${MOTIF_LENGTH}" \
        -p "${THREADS}" \
        -cache 2000
    
    log "Motif analysis complete for ${name}. Results in ${output_dir}/${name}"
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
        -i|--intersections)
            INTERSECTION_DIR="$2"
            shift 2
            ;;
        -s|--size)
            SIZE="$2"
            shift 2
            ;;
        -m|--motif-len)
            MOTIF_LENGTH="$2"
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

# Create output directory
create_directories

# Find all BED files in the intersection directory
bed_files=("${INTERSECTION_DIR}"/*.bed)

if [ ${#bed_files[@]} -eq 0 ]; then
    log "ERROR: No BED files found in ${INTERSECTION_DIR}"
    exit 1
fi

# Run motif analysis on each BED file
for bed_file in "${bed_files[@]}"; do
    run_motif_analysis "${bed_file}" "${OUTPUT_DIR}"
done

log "All motif analyses completed successfully."
echo ""
echo "Summary of analyses:"
echo "==================="
for bed_file in "${bed_files[@]}"; do
    name=$(basename "${bed_file}" .bed)
    echo "${name}: ${OUTPUT_DIR}/${name}"
done
echo ""

exit 0
