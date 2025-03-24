#!/usr/bin/env bash
#
# Gli3-Sufu Master Analysis Workflow
# 
# This script orchestrates a complete workflow for analyzing
# ChIP-seq, CUT&RUN, ATAC-seq, and Split-DamID data for
# Gli3-Sufu gene regulatory network analysis.
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
GENOME_DIR="${HOME}/Databases/${GENOME}"
OUTPUT_DIR="./gli3_sufu_results"
CHIPSEQ_SAMPLES=""
ATACSEQ_SAMPLES=""
CUTRUN_SAMPLES=""
SPLITDAMID_SAMPLES=""
THREADS=16

# Flags for which analyses to run
RUN_CHIPSEQ=false
RUN_ATACSEQ=false
RUN_CUTRUN=false
RUN_SPLITDAMID=false
RUN_INTEGRATION=true
RUN_DIFFBIND=true

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -o, --output DIR              Output directory (default: ${OUTPUT_DIR})"
    echo "  -g, --genome NAME             Genome name (default: ${GENOME})"
    echo "  --genome-dir DIR              Genome reference directory (default: ${GENOME_DIR})"
    echo "  -t, --threads NUM             Number of threads (default: ${THREADS})"
    echo ""
    echo "  --chipseq-samples FILE        Sample sheet for ChIP-seq analysis"
    echo "  --atacseq-samples FILE        Sample sheet for ATAC-seq analysis"
    echo "  --cutrun-samples FILE         Sample sheet for CUT&RUN analysis"
    echo "  --splitdamid-samples FILE     Sample sheet for Split-DamID analysis"
    echo ""
    echo "  --skip-chipseq                Skip ChIP-seq analysis"
    echo "  --skip-atacseq                Skip ATAC-seq analysis"
    echo "  --skip-cutrun                 Skip CUT&RUN analysis"
    echo "  --skip-splitdamid             Skip Split-DamID analysis"
    echo "  --skip-integration            Skip peak integration analysis"
    echo "  --skip-diffbind               Skip differential binding analysis"
    echo ""
    echo "  --only-chipseq                Run only ChIP-seq analysis"
    echo "  --only-atacseq                Run only ATAC-seq analysis"
    echo "  --only-cutrun                 Run only CUT&RUN analysis"
    echo "  --only-splitdamid             Run only Split-DamID analysis"
    echo "  --only-integration            Run only peak integration analysis"
    echo "  --only-diffbind               Run only differential binding analysis"
    echo ""
    echo "  -h, --help                    Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    local dependencies=(nextflow conda Rscript python3)
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
    
    log "All dependencies found."
}

check_genome_files() {
    log "Checking genome reference files..."
    
    # Required genome files
    local required_files=(
        "${GENOME_DIR}/${GENOME}.fa"
        "${GENOME_DIR}/${GENOME}.ncbiRefSeq.gtf"
        "${GENOME_DIR}/${GENOME}-blacklist.v2.bed"
    )
    
    local missing=()
    for file in "${required_files[@]}"; do
        if [ ! -f "${file}" ]; then
            missing+=("${file}")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        log "ERROR: Missing genome files: ${missing[*]}"
        log "Please ensure all required genome files are available."
        exit 1
    fi
    
    log "All required genome files found."
}

create_directories() {
    log "Creating output directories..."
    
    mkdir -p "${OUTPUT_DIR}"
    
    if [ "${RUN_CHIPSEQ}" = true ]; then
        mkdir -p "${OUTPUT_DIR}/chipseq"
    fi
    
    if [ "${RUN_ATACSEQ}" = true ]; then
        mkdir -p "${OUTPUT_DIR}/atacseq"
    fi
    
    if [ "${RUN_CUTRUN}" = true ]; then
        mkdir -p "${OUTPUT_DIR}/cutrun"
    fi
    
    if [ "${RUN_SPLITDAMID}" = true ]; then
        mkdir -p "${OUTPUT_DIR}/splitdamid"
    fi
    
    if [ "${RUN_INTEGRATION}" = true ]; then
        mkdir -p "${OUTPUT_DIR}/integration"
    fi
    
    if [ "${RUN_DIFFBIND}" = true ]; then
        mkdir -p "${OUTPUT_DIR}/diffbind"
    fi
    
    log "Output directories created."
}

run_chipseq_analysis() {
    log "Running ChIP-seq analysis..."
    
    # Check if sample sheet exists
    if [ ! -f "${CHIPSEQ_SAMPLES}" ]; then
        log "ERROR: ChIP-seq sample sheet not found: ${CHIPSEQ_SAMPLES}"
        return 1
    fi
    
    # Run nf-core chipseq pipeline
    log "Launching nf-core/chipseq workflow..."
    nextflow run nf-core/chipseq \
        -profile singularity \
        --input "${CHIPSEQ_SAMPLES}" \
        --genome "${GENOME}" \
        --fasta "${GENOME_DIR}/${GENOME}.fa" \
        --gtf "${GENOME_DIR}/${GENOME}.ncbiRefSeq.gtf" \
        --blacklist "${GENOME_DIR}/${GENOME}-blacklist.v2.bed" \
        --max_cpus "${THREADS}" \
        --outdir "${OUTPUT_DIR}/chipseq"
    
    log "ChIP-seq analysis completed."
}

run_atacseq_analysis() {
    log "Running ATAC-seq analysis..."
    
    # Check if sample sheet exists
    if [ ! -f "${ATACSEQ_SAMPLES}" ]; then
        log "ERROR: ATAC-seq sample sheet not found: ${ATACSEQ_SAMPLES}"
        return 1
    fi
    
    # Run ATAC-seq analysis with custom script
    log "Launching ATAC-seq workflow..."
    ./run-atacseq-workflow.sh \
        --sample-sheet "${ATACSEQ_SAMPLES}" \
        --output "${OUTPUT_DIR}/atacseq" \
        --genome "${GENOME}" \
        --fasta "${GENOME_DIR}/${GENOME}.fa" \
        --gtf "${GENOME_DIR}/${GENOME}.ncbiRefSeq.gtf" \
        --blacklist "${GENOME_DIR}/${GENOME}-blacklist.v2.bed" \
        --threads "${THREADS}"
    
    log "ATAC-seq analysis completed."
}

run_cutrun_analysis() {
    log "Running CUT&RUN analysis..."
    
    # Check if sample sheet exists
    if [ ! -f "${CUTRUN_SAMPLES}" ]; then
        log "ERROR: CUT&RUN sample sheet not found: ${CUTRUN_SAMPLES}"
        return 1
    fi
    
    # Run nf-core cutandrun pipeline
    log "Launching nf-core/cutandrun workflow..."
    nextflow run nf-core/cutandrun \
        -profile singularity \
        --input "${CUTRUN_SAMPLES}" \
        --genome "${GENOME}" \
        --fasta "${GENOME_DIR}/${GENOME}.fa" \
        --gtf "${GENOME_DIR}/${GENOME}.ncbiRefSeq.gtf" \
        --blacklist "${GENOME_DIR}/${GENOME}-blacklist.v2.bed" \
        --peakcaller seacr \
        --seacr_stringent relaxed \
        --max_cpus "${THREADS}" \
        --outdir "${OUTPUT_DIR}/cutrun"
    
    log "CUT&RUN analysis completed."
}

run_splitdamid_analysis() {
    log "Running Split-DamID analysis..."
    
    # Check if sample file exists
    if [ ! -f "${SPLITDAMID_SAMPLES}" ]; then
        log "ERROR: Split-DamID sample file not found: ${SPLITDAMID_SAMPLES}"
        return 1
    fi
    
    # Check if GATC sites file exists, create if not
    if [ ! -f "${GENOME_DIR}/${GENOME}_GATC_sites.bed" ]; then
        log "GATC sites file not found, creating it..."
        python3 create_gatc_sites.py "${GENOME_DIR}/${GENOME}.fa" "${GENOME_DIR}/${GENOME}_GATC_sites.bed"
    fi
    
    # Check if blacklist+mtDNA file exists, create if not
    if [ ! -f "${GENOME_DIR}/blacklist_mtDNA.bed" ]; then
        log "Blacklist+mtDNA file not found, creating it..."
        ./create_blacklist_mtdna.sh --genome "${GENOME}" --output "${GENOME_DIR}/blacklist_mtDNA.bed"
    fi
    
    # Run Split-DamID analysis with custom script
    log "Launching Split-DamID workflow..."
    ./run-splitdamid-workflow.sh \
        --samples "${SPLITDAMID_SAMPLES}" \
        --output "${OUTPUT_DIR}/splitdamid" \
        --genome "${GENOME}" \
        --fasta "${GENOME_DIR}/${GENOME}.fa" \
        --bowtie2 "${GENOME_DIR}/${GENOME}_bt2index" \
        --gatc "${GENOME_DIR}/${GENOME}_GATC_sites.bed" \
        --blacklist "${GENOME_DIR}/blacklist_mtDNA.bed" \
        --extension-script "extend_reads_to_GATC_revised.py" \
        --threads "${THREADS}"
    
    log "Split-DamID analysis completed."
}

run_peak_integration() {
    log "Running peak integration analysis..."
    
    # Prepare directories for peak files
    CHIP_PEAKS_DIR="${OUTPUT_DIR}/integration/chipseq_peaks"
    ATACSEQ_PEAKS_DIR="${OUTPUT_DIR}/integration/atacseq_peaks"
    CUTRUN_PEAKS_DIR="${OUTPUT_DIR}/integration/cutrun_peaks"
    SPLITDAMID_PEAKS_DIR="${OUTPUT_DIR}/integration/splitdamid_peaks"
    
    mkdir -p "${CHIP_PEAKS_DIR}" "${ATACSEQ_PEAKS_DIR}" "${CUTRUN_PEAKS_DIR}" "${SPLITDAMID_PEAKS_DIR}"
    
    # Copy and prepare ChIP-seq peak files
    if [ -d "${OUTPUT_DIR}/chipseq" ]; then
        log "Preparing ChIP-seq peak files for integration..."
        find "${OUTPUT_DIR}/chipseq" -name "*.narrowPeak" -o -name "*.broadPeak" -o -name "*consensus_peaks*.bed" | \
            while read peak_file; do
                base_name=$(basename "${peak_file}" | sed 's/\..*//')
                cp "${peak_file}" "${CHIP_PEAKS_DIR}/${base_name}.bed"
            done
    fi
    
    # Copy and prepare ATAC-seq peak files
    if [ -d "${OUTPUT_DIR}/atacseq" ]; then
        log "Preparing ATAC-seq peak files for integration..."
        find "${OUTPUT_DIR}/atacseq" -name "*consensus_peaks*.bed" | \
            while read peak_file; do
                base_name=$(basename "${peak_file}" | sed 's/\..*//')
                cp "${peak_file}" "${ATACSEQ_PEAKS_DIR}/${base_name}.bed"
            done
    fi
    
    # Copy and prepare CUT&RUN peak files
    if [ -d "${OUTPUT_DIR}/cutrun" ]; then
        log "Preparing CUT&RUN peak files for integration..."
        find "${OUTPUT_DIR}/cutrun" -name "*.narrowPeak" -o -name "*.broadPeak" -o -name "*consensus_peaks*.bed" | \
            while read peak_file; do
                base_name=$(basename "${peak_file}" | sed 's/\..*//')
                cp "${peak_file}" "${CUTRUN_PEAKS_DIR}/${base_name}.bed"
            done
    fi
    
    # Copy and prepare Split-DamID peak files
    if [ -d "${OUTPUT_DIR}/splitdamid" ]; then
        log "Preparing Split-DamID peak files for integration..."
        find "${OUTPUT_DIR}/splitdamid/peaks" -name "*broad_peaks.broadPeak" | \
            while read peak_file; do
                base_name=$(basename "${peak_file}" .broadPeak)
                cp "${peak_file}" "${SPLITDAMID_PEAKS_DIR}/${base_name}.bed"
            done
    fi
    
    # Run peak intersection analysis
    log "Running peak intersection workflow..."
    ./peak-intersection-workflow.sh \
        --output "${OUTPUT_DIR}/integration" \
        --genome "${GENOME}" \
        --chipseq-peaks "${CHIP_PEAKS_DIR}" \
        --atacseq-peaks "${ATACSEQ_PEAKS_DIR}" \
        --cutrun-peaks "${CUTRUN_PEAKS_DIR}" \
        --splitdamid-peaks "${SPLITDAMID_PEAKS_DIR}"
    
    log "Peak integration analysis completed."
}

run_diffbind_analysis() {
    log "Running differential binding analysis..."
    
    # Run DiffBind for each data type
    
    # ChIP-seq DiffBind
    if [ -d "${OUTPUT_DIR}/chipseq" ]; then
        log "Creating ChIP-seq DiffBind sample sheet..."
        ./create_diffbind_samplesheet.sh \
            --input "${OUTPUT_DIR}/chipseq" \
            --output "${OUTPUT_DIR}/diffbind/chipseq_samples.csv" \
            --tissue "ChIPseq"
        
        log "Running ChIP-seq DiffBind analysis..."
        Rscript splitdamid_diffbind_analysis.R \
            "${OUTPUT_DIR}/diffbind/chipseq_samples.csv" \
            "${OUTPUT_DIR}/diffbind/chipseq"
    fi
    
    # CUT&RUN DiffBind
    if [ -d "${OUTPUT_DIR}/cutrun" ]; then
        log "Creating CUT&RUN DiffBind sample sheet..."
        ./create_diffbind_samplesheet.sh \
            --input "${OUTPUT_DIR}/cutrun" \
            --output "${OUTPUT_DIR}/diffbind/cutrun_samples.csv" \
            --tissue "CUT&RUN"
        
        log "Running CUT&RUN DiffBind analysis..."
        Rscript splitdamid_diffbind_analysis.R \
            "${OUTPUT_DIR}/diffbind/cutrun_samples.csv" \
            "${OUTPUT_DIR}/diffbind/cutrun"
    fi
    
    # Split-DamID DiffBind
    if [ -d "${OUTPUT_DIR}/splitdamid" ]; then
        log "Creating Split-DamID DiffBind sample sheet..."
        ./create_diffbind_samplesheet.sh \
            --input "${OUTPUT_DIR}/splitdamid" \
            --output "${OUTPUT_DIR}/diffbind/splitdamid_samples.csv" \
            --tissue "NIH3T3"
        
        log "Running Split-DamID DiffBind analysis..."
        Rscript splitdamid_diffbind_analysis.R \
            "${OUTPUT_DIR}/diffbind/splitdamid_samples.csv" \
            "${OUTPUT_DIR}/diffbind/splitdamid"
    fi
    
    log "Differential binding analysis completed."
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
        --genome-dir)
            GENOME_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --chipseq-samples)
            CHIPSEQ_SAMPLES="$2"
            RUN_CHIPSEQ=true
            shift 2
            ;;
        --atacseq-samples)
            ATACSEQ_SAMPLES="$2"
            RUN_ATACSEQ=true
            shift 2
            ;;
        --cutrun-samples)
            CUTRUN_SAMPLES="$2"
            RUN_CUTRUN=true
            shift 2
            ;;
        --splitdamid-samples)
            SPLITDAMID_SAMPLES="$2"
            RUN_SPLITDAMID=true
            shift 2
            ;;
        --skip-chipseq)
            RUN_CHIPSEQ=false
            shift
            ;;
        --skip-atacseq)
            RUN_ATACSEQ=false
            shift
            ;;
        --skip-cutrun)
            RUN_CUTRUN=false
            shift
            ;;
        --skip-splitdamid)
            RUN_SPLITDAMID=false
            shift
            ;;
        --skip-integration)
            RUN_INTEGRATION=false
            shift
            ;;
        --skip-diffbind)
            RUN_DIFFBIND=false
            shift
            ;;
        --only-chipseq)
            RUN_CHIPSEQ=true
            RUN_ATACSEQ=false
            RUN_CUTRUN=false
            RUN_SPLITDAMID=false
            RUN_INTEGRATION=false
            RUN_DIFFBIND=false
            shift
            ;;
        --only-atacseq)
            RUN_CHIPSEQ=false
            RUN_ATACSEQ=true
            RUN_CUTRUN=false
            RUN_SPLITDAMID=false
            RUN_INTEGRATION=false
            RUN_DIFFBIND=false
            shift
            ;;
        --only-cutrun)
            RUN_CHIPSEQ=false
            RUN_ATACSEQ=false
            RUN_CUTRUN=true
            RUN_SPLITDAMID=false
            RUN_INTEGRATION=false
            RUN_DIFFBIND=false
            shift
            ;;
        --only-splitdamid)
            RUN_CHIPSEQ=false
            RUN_ATACSEQ=false
            RUN_CUTRUN=false
            RUN_SPLITDAMID=true
            RUN_INTEGRATION=false
            RUN_DIFFBIND=false
            shift
            ;;
        --only-integration)
            RUN_CHIPSEQ=false
            RUN_ATACSEQ=false
            RUN_CUTRUN=false
            RUN_SPLITDAMID=false
            RUN_INTEGRATION=true
            RUN_DIFFBIND=false
            shift
            ;;
        --only-diffbind)
            RUN_CHIPSEQ=false
            RUN_ATACSEQ=false
            RUN_CUTRUN=false
            RUN_SPLITDAMID=false
            RUN_INTEGRATION=false
            RUN_DIFFBIND=true
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

# Print analysis configuration
log "Gli3-Sufu Analysis Workflow Configuration:"
log "  Genome: ${GENOME}"
log "  Genome Directory: ${GENOME_DIR}"
log "  Output Directory: ${OUTPUT_DIR}"
log "  Threads: ${THREADS}"
log "  Analyses:"
log "    - ChIP-seq: $(if [ "${RUN_CHIPSEQ}" = true ]; then echo "YES"; else echo "NO"; fi)"
log "    - ATAC-seq: $(if [ "${RUN_ATACSEQ}" = true ]; then echo "YES"; else echo "NO"; fi)"
log "    - CUT&RUN: $(if [ "${RUN_CUTRUN}" = true ]; then echo "YES"; else echo "NO"; fi)"
log "    - Split-DamID: $(if [ "${RUN_SPLITDAMID}" = true ]; then echo "YES"; else echo "NO"; fi)"
log "    - Peak Integration: $(if [ "${RUN_INTEGRATION}" = true ]; then echo "YES"; else echo "NO"; fi)"
log "    - Differential Binding: $(if [ "${RUN_DIFFBIND}" = true ]; then echo "YES"; else echo "NO"; fi)"

# Check dependencies
check_dependencies

# Check genome files
check_genome_files

# Create output directories
create_directories

# Run selected analyses
if [ "${RUN_CHIPSEQ}" = true ]; then
    run_chipseq_analysis
fi

if [ "${RUN_ATACSEQ}" = true ]; then
    run_atacseq_analysis
fi

if [ "${RUN_CUTRUN}" = true ]; then
    run_cutrun_analysis
fi

if [ "${RUN_SPLITDAMID}" = true ]; then
    run_splitdamid_analysis
fi

if [ "${RUN_INTEGRATION}" = true ]; then
    run_peak_integration
fi

if [ "${RUN_DIFFBIND}" = true ]; then
    run_diffbind_analysis
fi

log "Gli3-Sufu analysis workflow completed successfully!"
log "Results can be found in: ${OUTPUT_DIR}"

exit 0