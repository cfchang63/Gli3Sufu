#!/usr/bin/env bash
#
# Create a combined blacklist and mtDNA exclusion file for Split-DamID analysis
# 
# This script downloads the ENCODE blacklist for mm10 and adds the mitochondrial
# chromosome to create a combined exclusion list.
#
# Author: Ching-Fang Chang
# Date: March 23, 2025
#

set -e  # Exit on error

# CONFIGURATION
##############################################################################

# Default parameters (can be overridden via command line arguments)
GENOME="mm10"
OUTPUT_FILE="blacklist_mtDNA.bed"

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -g, --genome NAME        Genome name (mm10 or hg38, default: ${GENOME})"
    echo "  -o, --output FILE        Output file name (default: ${OUTPUT_FILE})"
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
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
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

# Download appropriate blacklist
log "Downloading ${GENOME} blacklist from ENCODE..."

if [ "${GENOME}" = "mm10" ]; then
    # Mouse mm10 blacklist
    wget -O ${GENOME}-blacklist.v2.bed.gz https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz
elif [ "${GENOME}" = "hg38" ]; then
    # Human hg38 blacklist
    wget -O ${GENOME}-blacklist.v2.bed.gz https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
else
    log "ERROR: Unsupported genome: ${GENOME}. Only mm10 and hg38 are supported."
    exit 1
fi

# Uncompress blacklist
log "Uncompressing blacklist..."
gunzip -f ${GENOME}-blacklist.v2.bed.gz

# Create mitochondrial exclusion region
log "Creating mitochondrial exclusion region..."

if [ "${GENOME}" = "mm10" ]; then
    # Mouse mitochondrial chromosome length is 16,299 bp
    echo -e "chrM\t0\t16299\tMitochondrial" > mtDNA_regions.bed
elif [ "${GENOME}" = "hg38" ]; then
    # Human mitochondrial chromosome length is 16,569 bp
    echo -e "chrM\t0\t16569\tMitochondrial" > mtDNA_regions.bed
fi

# Combine blacklist and mtDNA
log "Combining blacklist and mtDNA exclusion regions..."
cat ${GENOME}-blacklist.v2.bed mtDNA_regions.bed > ${OUTPUT_FILE}

# Count regions
TOTAL_REGIONS=$(wc -l < ${OUTPUT_FILE})
log "Created combined exclusion file with ${TOTAL_REGIONS} regions"

# Clean up
log "Cleaning up temporary files..."
rm mtDNA_regions.bed

log "Success! Combined blacklist and mtDNA exclusion file created: ${OUTPUT_FILE}"
exit 0
