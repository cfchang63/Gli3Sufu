#!/usr/bin/env bash
#
# ChIP-seq Analysis Workflow
# 
# This script automates the process of running ChIP-seq analysis using the 
# nf-core/chipseq pipeline and performs post-processing for peak integration
# with other data types.
#
# Author: Ching-Fang Chang
# Date: March 23, 2025
#

Usage:
## Basic usage with default parameters
# ./chipseq-workflow.sh --sample-sheet samples.csv --output chipseq_results

## Single-end data with narrow peaks
# ./chipseq-workflow.sh --sample-sheet samples.csv --output chipseq_results --single-end --read-length 75

## Integrate with external peaks (e.g., CUT&RUN peaks)
# ./chipseq-workflow.sh --sample-sheet samples.csv --output chipseq_results --integration-peaks gli3_peaks.bed

set -e  # Exit on error
set -o pipefail  # Exit if any command in pipe fails

# CONFIGURATION
##############################################################################

# Default parameters (can be overridden via command line arguments)
GENOME="mm10"
GENOME_DIR="${HOME}/Databases/${GENOME}"
OUTPUT_DIR="./chipseq_results"
FASTA="${GENOME_DIR}/${GENOME}.fa"
GTF="${GENOME_DIR}/${GENOME}.ncbiRefSeq.gtf"
BLACKLIST="${GENOME_DIR}/${GENOME}-blacklist.v2.bed"
BOWTIE2_INDEX="${GENOME_DIR}"
CHROMAP_INDEX="${GENOME_DIR}/${GENOME}.chrom.sizes"
MACS_GSIZE="2652783500.0"  # mm10 genome size
SAMPLE_SHEET=""
THREADS=16
SINGLE_END=false
READ_LENGTH="75"
NARROW_PEAK=true
INTEGRATION_PEAKS=""

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -s, --sample-sheet FILE  Sample sheet for nf-core/chipseq (required)"
    echo "  -o, --output DIR         Output directory (default: ${OUTPUT_DIR})"
    echo "  -g, --genome NAME        Genome name (default: ${GENOME})"
    echo "  --genome-dir DIR         Genome reference directory (default: ${GENOME_DIR})"
    echo "  --fasta FILE             Genome FASTA file (default: ${FASTA})"
    echo "  --gtf FILE               GTF file path (default: ${GTF})"
    echo "  --blacklist FILE         Blacklist regions BED file (default: ${BLACKLIST})"
    echo "  --bowtie2-index DIR      Bowtie2 index directory (default: ${BOWTIE2_INDEX})"
    echo "  --chromap-index FILE     Chromosome sizes file (default: ${CHROMAP_INDEX})"
    echo "  --macs-gsize SIZE        Genome size for MACS2 (default: ${MACS_GSIZE})"
    echo "  -t, --threads NUM        Number of threads (default: ${THREADS})"
    echo "  --single-end             Specify for single-end sequencing data (default: paired-end)"
    echo "  --read-length NUM        Read length (default: ${READ_LENGTH})"
    echo "  --broad-peak             Call broad peaks instead of narrow peaks"
    echo "  --integration-peaks FILE Path to BED file for peak integration (optional)"
    echo "  -h, --help               Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    local dependencies=(conda)
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
    
    # Check if conda environment nf-core exists
    if ! conda env list | grep -q "nf-core"; then
        log "Conda environment 'nf-core' not found. Creating it..."
        conda create --yes --name nf-core python=3.12 nf-core nextflow
    else
        log "Found conda environment 'nf-core'"
    fi
    
    log "All dependencies found."
}

create_directories() {
    log "Creating output directories..."
    
    mkdir -p "${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}/integration"
    mkdir -p "${OUTPUT_DIR}/annotation"
    
    log "Output directories created."
}

check_sample_sheet() {
    log "Checking sample sheet: ${SAMPLE_SHEET}"
    
    if [ ! -f "${SAMPLE_SHEET}" ]; then
        log "ERROR: Sample sheet not found: ${SAMPLE_SHEET}"
        exit 1
    fi
    
    # Check if the sample sheet has the correct format
    local header=$(head -n 1 "${SAMPLE_SHEET}")
    if [[ ! $header == "sample,fastq_1,fastq_2,antibody,control" ]]; then
        log "ERROR: Invalid sample sheet format. Header should be 'sample,fastq_1,fastq_2,antibody,control'"
        exit 1
    fi
    
    # Check if all fastq files exist
    local fastq_cols="2,3"
    if [ "${SINGLE_END}" = true ]; then
        fastq_cols="2"
    fi
    
    local missing_files=()
    while IFS=, read -r sample fastq1 fastq2 antibody control; do
        # Skip header
        if [ "${sample}" = "sample" ]; then
            continue
        fi
        
        # Check fastq_1
        if [ ! -f "${fastq1}" ] && [ -n "${fastq1}" ]; then
            missing_files+=("${fastq1}")
        fi
        
        # Check fastq_2 if paired-end
        if [ "${SINGLE_END}" = false ] && [ ! -f "${fastq2}" ] && [ -n "${fastq2}" ]; then
            missing_files+=("${fastq2}")
        fi
    done < "${SAMPLE_SHEET}"
    
    if [ ${#missing_files[@]} -gt 0 ]; then
        log "WARNING: Some FASTQ files not found:"
        for file in "${missing_files[@]}"; do
            log "  ${file}"
        done
        read -p "Continue anyway? (y/n) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
    fi
    
    log "Sample sheet format is valid."
}

run_nf_core_chipseq() {
    log "Running nf-core/chipseq pipeline..."
    
    # Activate nf-core conda environment
    log "Activating conda environment 'nf-core'..."
    eval "$(conda shell.bash hook)"
    conda activate nf-core
    
    # Prepare command
    local cmd="nextflow run nf-core/chipseq"
    cmd+=" --input ${SAMPLE_SHEET}"
    cmd+=" --outdir ${OUTPUT_DIR}"
    cmd+=" --aligner bowtie2"
    cmd+=" --fasta ${FASTA}"
    cmd+=" --gtf ${GTF}"
    cmd+=" --blacklist ${BLACKLIST}"
    cmd+=" --bowtie2_index ${BOWTIE2_INDEX}"
    cmd+=" --chromap_index ${CHROMAP_INDEX}"
    cmd+=" --macs_gsize ${MACS_GSIZE}"
    cmd+=" --read_length ${READ_LENGTH}"
    cmd+=" -profile singularity"
    cmd+=" -with-report ${OUTPUT_DIR}/chipseq_output_report.html"
    cmd+=" -with-timeline ${OUTPUT_DIR}/chipseq_output_timeline.html"
    cmd+=" --params-out ${OUTPUT_DIR}/params-output.json"
    cmd+=" --max_cpus ${THREADS}"
    
    # Add single-end flag if needed
    if [ "${SINGLE_END}" = true ]; then
        cmd+=" --single_end"
    fi
    
    # Add narrow/broad peak flag
    if [ "${NARROW_PEAK}" = true ]; then
        cmd+=" --narrow_peak"
    else
        cmd+=" --broad_peak"
    fi
    
    # Run the pipeline
    log "Executing command: ${cmd}"
    eval "${cmd}"
    
    # Deactivate the environment
    conda deactivate
    
    log "nf-core/chipseq pipeline completed."
}

extract_antibodies() {
    log "Extracting antibodies from the sample sheet..."
    
    local antibodies=()
    
    # Skip header line and extract unique antibodies
    antibodies=$(tail -n +2 "${SAMPLE_SHEET}" | cut -d',' -f4 | sort | uniq | grep -v '^$')
    
    log "Found antibodies: ${antibodies}"
    echo "${antibodies}"
}

integrate_with_external_peaks() {
    local integration_peaks="$1"
    local antibody="$2"
    
    log "Integrating ${antibody} peaks with external peaks..."
    
    if [ ! -f "${integration_peaks}" ]; then
        log "ERROR: Integration peaks file not found: ${integration_peaks}"
        return 1
    fi
    
    # Determine the path to consensus peaks based on narrow/broad peak setting
    local peak_type="narrowPeak"
    if [ "${NARROW_PEAK}" = false ]; then
        peak_type="broadPeak"
    fi
    
    local consensus_peaks="${OUTPUT_DIR}/bowtie2/mergedLibrary/macs2/${peak_type}/consensus/${antibody}/${antibody}.consensus_peaks.bed"
    
    if [ ! -f "${consensus_peaks}" ]; then
        log "WARNING: Consensus peaks file not found at ${consensus_peaks}"
        log "Searching for alternative location..."
        
        # Try to find consensus peaks file in alternative locations
        local alternative_peaks=$(find "${OUTPUT_DIR}" -name "${antibody}.consensus_peaks.bed" | head -n 1)
        
        if [ -z "${alternative_peaks}" ]; then
            log "ERROR: Could not find consensus peaks file for ${antibody}"
            return 1
        else
            consensus_peaks="${alternative_peaks}"
            log "Found consensus peaks at ${consensus_peaks}"
        fi
    fi
    
    # Prepare sorted BED files
    log "Preparing sorted BED files..."
    
    # Extract columns 1-3 from consensus peaks (chromosome, start, end)
    cut -f1-3 "${consensus_peaks}" > "${OUTPUT_DIR}/integration/${antibody}_peaks.bed"
    sort -k1,1 -k2,2n "${OUTPUT_DIR}/integration/${antibody}_peaks.bed" > "${OUTPUT_DIR}/integration/${antibody}_peaks_sorted.bed"
    
    # Sort external peaks
    sort -k1,1 -k2,2n "${integration_peaks}" > "${OUTPUT_DIR}/integration/external_peaks_sorted.bed"
    
    # Find overlapping peaks
    log "Finding overlapping peaks..."
    bedtools intersect -a "${OUTPUT_DIR}/integration/${antibody}_peaks_sorted.bed" \
        -b "${OUTPUT_DIR}/integration/external_peaks_sorted.bed" \
        -wa -wb > "${OUTPUT_DIR}/integration/${antibody}_overlapping_peaks.bed"
    
    # Annotate overlapping peaks with HOMER
    log "Annotating overlapping peaks with HOMER..."
    annotatePeaks.pl "${OUTPUT_DIR}/integration/${antibody}_overlapping_peaks.bed" \
        "${GENOME}" > "${OUTPUT_DIR}/annotation/${antibody}_overlapping_peaks_annotation.txt"
    
    log "Integration and annotation complete for ${antibody}"
    log "Results saved to ${OUTPUT_DIR}/annotation/${antibody}_overlapping_peaks_annotation.txt"
}

process_peak_results() {
    log "Processing peak results..."
    
    # Get list of antibodies
    local antibodies=$(extract_antibodies)
    
    for antibody in ${antibodies}; do
        if [ -z "${antibody}" ] || [ "${antibody}" = "Input" ]; then
            continue
        fi
        
        log "Processing results for antibody: ${antibody}"
        
        # Determine the path to consensus peaks and annotation
        local peak_type="narrowPeak"
        if [ "${NARROW_PEAK}" = false ]; then
            peak_type="broadPeak"
        fi
        
        local consensus_peaks="${OUTPUT_DIR}/bowtie2/mergedLibrary/macs2/${peak_type}/consensus/${antibody}/${antibody}.consensus_peaks.bed"
        local annotation_file="${OUTPUT_DIR}/bowtie2/mergedLibrary/macs2/${peak_type}/consensus/${antibody}/${antibody}.consensus_peaks.annotatePeaks.txt"
        
        # Copy annotation file to annotation directory
        if [ -f "${annotation_file}" ]; then
            cp "${annotation_file}" "${OUTPUT_DIR}/annotation/${antibody}_annotation.txt"
            log "Copied annotation file to ${OUTPUT_DIR}/annotation/${antibody}_annotation.txt"
        else
            log "WARNING: Annotation file not found at ${annotation_file}"
        }
        
        # If integration peaks are provided, integrate them
        if [ -n "${INTEGRATION_PEAKS}" ]; then
            integrate_with_external_peaks "${INTEGRATION_PEAKS}" "${antibody}"
        fi
    done
    
    log "Peak processing completed."
}

create_sample_sheet_template() {
    local template_file="chipseq_sample_template.csv"
    
    log "Creating sample sheet template at ${template_file}..."
    
    cat > "${template_file}" << EOF
sample,fastq_1,fastq_2,antibody,control
WT_Input,/path/to/WT_Input_R1.fastq.gz,/path/to/WT_Input_R2.fastq.gz,,
WT_TF1_Rep1,/path/to/WT_TF1_Rep1_R1.fastq.gz,/path/to/WT_TF1_Rep1_R2.fastq.gz,TF1,WT_Input
WT_TF1_Rep2,/path/to/WT_TF1_Rep2_R1.fastq.gz,/path/to/WT_TF1_Rep2_R2.fastq.gz,TF1,WT_Input
WT_TF2_Rep1,/path/to/WT_TF2_Rep1_R1.fastq.gz,/path/to/WT_TF2_Rep1_R2.fastq.gz,TF2,WT_Input
WT_TF2_Rep2,/path/to/WT_TF2_Rep2_R1.fastq.gz,/path/to/WT_TF2_Rep2_R2.fastq.gz,TF2,WT_Input
EOF
    
    if [ "${SINGLE_END}" = true ]; then
        log "Note: For single-end data, leave the fastq_2 column empty"
        cat > "${template_file}" << EOF
sample,fastq_1,fastq_2,antibody,control
WT_Input,/path/to/WT_Input.fastq.gz,,Input,
WT_TF1_Rep1,/path/to/WT_TF1_Rep1.fastq.gz,,TF1,WT_Input
WT_TF1_Rep2,/path/to/WT_TF1_Rep2.fastq.gz,,TF1,WT_Input
WT_TF2_Rep1,/path/to/WT_TF2_Rep1.fastq.gz,,TF2,WT_Input
WT_TF2_Rep2,/path/to/WT_TF2_Rep2.fastq.gz,,TF2,WT_Input
EOF
    }
    
    log "Sample sheet template created. Edit this file with your sample information."
    log "Then run the script again with -s ${template_file}"
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
        --genome-dir)
            GENOME_DIR="$2"
            shift 2
            ;;
        --fasta)
            FASTA="$2"
            shift 2
            ;;
        --gtf)
            GTF="$2"
            shift 2
            ;;
        --blacklist)
            BLACKLIST="$2"
            shift 2
            ;;
        --bowtie2-index)
            BOWTIE2_INDEX="$2"
            shift 2
            ;;
        --chromap-index)
            CHROMAP_INDEX="$2"
            shift 2
            ;;
        --macs-gsize)
            MACS_GSIZE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --single-end)
            SINGLE_END=true
            shift
            ;;
        --read-length)
            READ_LENGTH="$2"
            shift 2
            ;;
        --broad-peak)
            NARROW_PEAK=false
            shift
            ;;
        --integration-peaks)
            INTEGRATION_PEAKS="$2"
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

# If no sample sheet is provided, create a template and exit
if [ -z "${SAMPLE_SHEET}" ]; then
    log "No sample sheet provided."
    create_sample_sheet_template
    exit 0
fi

# Check dependencies
check_dependencies

# Create output directories
create_directories

# Check sample sheet format and existence of FASTQ files
check_sample_sheet

# Run nf-core/chipseq pipeline
run_nf_core_chipseq

# Process peak results
process_peak_results

log "ChIP-seq analysis workflow completed successfully!"
log "Results can be found in: ${OUTPUT_DIR}"
log "Annotation files: ${OUTPUT_DIR}/annotation/"

# Print summary
echo ""
echo "Summary of results:"
echo "==================="

# Print antibodies processed
antibodies=$(extract_antibodies)
for antibody in ${antibodies}; do
    if [ -z "${antibody}" ] || [ "${antibody}" = "Input" ]; then
        continue
    fi
    
    anno_file="${OUTPUT_DIR}/annotation/${antibody}_annotation.txt"
    if [ -f "${anno_file}" ]; then
        peak_count=$(wc -l < "${anno_file}")
        peak_count=$((peak_count - 1))  # Subtract header line
        echo "${antibody}: ${peak_count} peaks"
    fi
    
    if [ -n "${INTEGRATION_PEAKS}" ]; then
        overlap_file="${OUTPUT_DIR}/annotation/${antibody}_overlapping_peaks_annotation.txt"
        if [ -f "${overlap_file}" ]; then
            overlap_count=$(wc -l < "${overlap_file}")
            overlap_count=$((overlap_count - 1))  # Subtract header line
            echo "${antibody} overlapping peaks: ${overlap_count}"
        fi
    fi
done

echo ""
log "Done."
exit 0
