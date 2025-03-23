#!/usr/bin/env bash
#
# ATAC-seq Analysis Pipeline
# 
# This script processes ATAC-seq data through a complete workflow:
# - Quality control
# - Adapter trimming
# - Alignment
# - Post-alignment filtering
# - Peak calling
# - Peak analysis
#
# Author: [Your Name]
# Date: March 23, 2025
# License: [License]
#

set -e  # Exit on error
set -o pipefail  # Exit if any command in pipe fails

# CONFIGURATION
##############################################################################

# Set default parameters (can be overridden via command line arguments)
THREADS=24
GENOME="mm10"
GENOME_INDEX="${HOME}/Databases/${GENOME}/${GENOME}_bt2index"
BLACKLIST="${HOME}/Databases/${GENOME}/${GENOME}-blacklist.v2.bed"
CHROM_SIZES="${HOME}/Databases/${GENOME}/${GENOME}.chrom.sizes"
OUTPUT_DIR="./results"
SAMPLE_FILE="samples.txt"

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -s, --samples FILE     Text file with sample names (default: ${SAMPLE_FILE})"
    echo "  -o, --output DIR       Output directory (default: ${OUTPUT_DIR})"
    echo "  -g, --genome NAME      Genome name (default: ${GENOME})"
    echo "  -t, --threads NUM      Number of threads (default: ${THREADS})"
    echo "  -h, --help             Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    local dependencies=(fastqc trimmomatic bowtie2 samtools bedtools picard alignmentSieve macs3 bedGraphToBigWig)
    local missing=()
    
    log "Checking dependencies..."
    for cmd in "${dependencies[@]}"; do
        if ! command -v ${cmd} &> /dev/null; then
            missing+=(${cmd})
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        log "ERROR: Missing dependencies: ${missing[*]}"
        log "Please install missing dependencies and try again."
        exit 1
    fi
    
    log "All dependencies found."
}

check_files_exist() {
    local file_list=("$@")
    local missing=()
    
    for file in "${file_list[@]}"; do
        if [ ! -f "${file}" ]; then
            missing+=(${file})
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        log "ERROR: Missing files: ${missing[*]}"
        exit 1
    fi
}

create_directories() {
    local directories=(
        "${OUTPUT_DIR}"
        "${OUTPUT_DIR}/fastqc"
        "${OUTPUT_DIR}/trimmed"
        "${OUTPUT_DIR}/aligned"
        "${OUTPUT_DIR}/filtered"
        "${OUTPUT_DIR}/peaks"
        "${OUTPUT_DIR}/bigwig"
        "${OUTPUT_DIR}/differential"
    )
    
    log "Creating output directories..."
    for dir in "${directories[@]}"; do
        mkdir -p "${dir}"
    done
}

rename_fastq_files() {
    log "Standardizing FASTQ filenames..."
    
    # Remove sequencer info from filenames
    for file in *.fastq.gz; do
        new_name=$(echo ${file} | sed 's/_S[0-9]*_R/_R/')
        if [ "${file}" != "${new_name}" ]; then
            mv "${file}" "${new_name}"
        fi
    done
    
    # Remove _001 suffix
    for file in *.fastq.gz; do
        new_name=$(echo ${file} | sed 's/_001\././')
        if [ "${file}" != "${new_name}" ]; then
            mv "${file}" "${new_name}"
        fi
    done
    
    log "Filename standardization complete."
}

run_fastqc() {
    local sample=$1
    
    log "Running FastQC on ${sample}..."
    fastqc "${sample}_R1.fastq.gz" -d . -o "${OUTPUT_DIR}/fastqc"
    fastqc "${sample}_R2.fastq.gz" -d . -o "${OUTPUT_DIR}/fastqc"
    log "FastQC complete for ${sample}."
}

run_trimming() {
    local sample=$1
    
    log "Trimming adapters for ${sample}..."
    trimmomatic PE -threads ${THREADS} \
        "${sample}_R1.fastq.gz" \
        "${sample}_R2.fastq.gz" \
        "${OUTPUT_DIR}/trimmed/trimmed_${sample}_R1.fastq.gz" \
        "${OUTPUT_DIR}/trimmed/unpaired_${sample}_R1.fastq.gz" \
        "${OUTPUT_DIR}/trimmed/trimmed_${sample}_R2.fastq.gz" \
        "${OUTPUT_DIR}/trimmed/unpaired_${sample}_R2.fastq.gz" \
        LEADING:5 TRAILING:5 MINLEN:50
    
    log "Trimming complete for ${sample}."
}

run_alignment() {
    local sample=$1
    
    log "Aligning ${sample} to ${GENOME}..."
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -X 2000 -p ${THREADS} \
        -x "${GENOME_INDEX}" \
        -1 "${OUTPUT_DIR}/trimmed/trimmed_${sample}_R1.fastq.gz" \
        -2 "${OUTPUT_DIR}/trimmed/trimmed_${sample}_R2.fastq.gz" | \
        samtools view -bS - > "${OUTPUT_DIR}/aligned/${sample}.bam"
    
    log "Sorting BAM file for ${sample}..."
    samtools sort "${OUTPUT_DIR}/aligned/${sample}.bam" -o "${OUTPUT_DIR}/aligned/${sample}_sorted.bam"
    
    log "Indexing BAM file for ${sample}..."
    samtools index "${OUTPUT_DIR}/aligned/${sample}_sorted.bam"
    
    log "Alignment complete for ${sample}."
}

filter_mitochondrial() {
    local sample=$1
    
    log "Generating alignment statistics for ${sample}..."
    samtools idxstats "${OUTPUT_DIR}/aligned/${sample}_sorted.bam" > "${OUTPUT_DIR}/aligned/${sample}_sorted.idxstats"
    samtools flagstat "${OUTPUT_DIR}/aligned/${sample}_sorted.bam" > "${OUTPUT_DIR}/aligned/${sample}_sorted.flagstat"
    
    log "Removing mitochondrial reads for ${sample}..."
    samtools view -h "${OUTPUT_DIR}/aligned/${sample}_sorted.bam" | \
        grep -v chrM | \
        samtools sort -O bam -o "${OUTPUT_DIR}/filtered/${sample}.rmChrM.bam"
    
    log "Mitochondrial filtering complete for ${sample}."
}

remove_duplicates() {
    local sample=$1
    
    log "Adding read groups to ${sample}..."
    picard AddOrReplaceReadGroups \
        I="${OUTPUT_DIR}/filtered/${sample}.rmChrM.bam" \
        O="${OUTPUT_DIR}/filtered/${sample}.withRG.bam" \
        RGID=id \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=sample
    
    log "Marking duplicates for ${sample}..."
    picard MarkDuplicates QUIET=true \
        INPUT="${OUTPUT_DIR}/filtered/${sample}.withRG.bam" \
        OUTPUT="${OUTPUT_DIR}/filtered/${sample}.marked.bam" \
        METRICS_FILE="${OUTPUT_DIR}/filtered/${sample}.dup.metrics" \
        REMOVE_DUPLICATES=false \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT
    
    log "Filtering reads for ${sample}..."
    samtools view -h -b -f 2 -F 1548 -q 30 "${OUTPUT_DIR}/filtered/${sample}.marked.bam" | \
        samtools sort -o "${OUTPUT_DIR}/filtered/${sample}.filtered.bam"
    
    samtools index "${OUTPUT_DIR}/filtered/${sample}.filtered.bam"
    
    log "Duplicate removal complete for ${sample}."
}

filter_blacklist() {
    local sample=$1
    
    log "Filtering blacklisted regions for ${sample}..."
    bedtools intersect -nonamecheck -v -abam "${OUTPUT_DIR}/filtered/${sample}.filtered.bam" \
        -b "${BLACKLIST}" > "${OUTPUT_DIR}/filtered/${sample}.tmp.bam"
    
    samtools sort -O bam -o "${OUTPUT_DIR}/filtered/${sample}.blacklist-filtered.bam" \
        "${OUTPUT_DIR}/filtered/${sample}.tmp.bam"
    
    samtools index "${OUTPUT_DIR}/filtered/${sample}.blacklist-filtered.bam"
    
    rm "${OUTPUT_DIR}/filtered/${sample}.tmp.bam"
    
    log "Blacklist filtering complete for ${sample}."
}

shift_reads() {
    local sample=$1
    
    log "Shifting read coordinates for ${sample}..."
    alignmentSieve --numberOfProcessors max --ATACshift \
        --blackListFileName "${BLACKLIST}" \
        --bam "${OUTPUT_DIR}/filtered/${sample}.blacklist-filtered.bam" \
        -o "${OUTPUT_DIR}/filtered/${sample}.shifted.bam"
    
    samtools sort -o "${OUTPUT_DIR}/filtered/${sample}.shifted.sorted.bam" \
        "${OUTPUT_DIR}/filtered/${sample}.shifted.bam"
    
    samtools index "${OUTPUT_DIR}/filtered/${sample}.shifted.sorted.bam"
    
    log "Read shifting complete for ${sample}."
}

call_peaks() {
    local sample=$1
    
    log "Converting BAM to BED for ${sample}..."
    macs3 randsample -i "${OUTPUT_DIR}/filtered/${sample}.shifted.sorted.bam" \
        -f BAMPE -p 100 -o "${OUTPUT_DIR}/peaks/${sample}.bed"
    
    log "Calling peaks for ${sample}..."
    macs3 callpeak -f BEDPE --nomodel --shift -37 --extsize 73 -g 2652783500 \
        -B --broad --keep-dup all --cutoff-analysis \
        -n "${sample}" \
        -t "${OUTPUT_DIR}/peaks/${sample}.bed" \
        --outdir "${OUTPUT_DIR}/peaks/${sample}" \
        2> "${OUTPUT_DIR}/peaks/${sample}.macs3.log"
    
    log "Peak calling complete for ${sample}."
}

generate_bigwigs() {
    local sample=$1
    
    log "Generating coverage files for ${sample}..."
    bamCoverage --bam "${OUTPUT_DIR}/filtered/${sample}.shifted.bam" \
        -o "${OUTPUT_DIR}/bigwig/${sample}.bw" \
        --binSize 10 \
        --normalizeUsing CPM \
        --blackListFileName "${BLACKLIST}"
    
    log "Coverage generation complete for ${sample}."
}

generate_enrichment() {
    local sample=$1
    
    log "Generating fold enrichment bedGraph for ${sample}..."
    macs3 bdgcmp \
        -t "${OUTPUT_DIR}/peaks/${sample}/${sample}_treat_pileup.bdg" \
        -c "${OUTPUT_DIR}/peaks/${sample}/${sample}_control_lambda.bdg" \
        -m FE \
        -o "${OUTPUT_DIR}/peaks/${sample}_FE.bdg"
    
    macs3 bdgcmp \
        -t "${OUTPUT_DIR}/peaks/${sample}/${sample}_treat_pileup.bdg" \
        -c "${OUTPUT_DIR}/peaks/${sample}/${sample}_control_lambda.bdg" \
        -m ppois \
        -o "${OUTPUT_DIR}/peaks/${sample}_ppois.bdg"
    
    log "Sorting bedGraph files for ${sample}..."
    sort -k1,1 -k2,2n "${OUTPUT_DIR}/peaks/${sample}_FE.bdg" > "${OUTPUT_DIR}/peaks/${sample}_FE.sorted.bdg"
    sort -k1,1 -k2,2n "${OUTPUT_DIR}/peaks/${sample}_ppois.bdg" > "${OUTPUT_DIR}/peaks/${sample}_ppois.sorted.bdg"
    
    log "Converting bedGraph to bigWig for ${sample}..."
    bedGraphToBigWig "${OUTPUT_DIR}/peaks/${sample}_FE.sorted.bdg" \
        "${CHROM_SIZES}" \
        "${OUTPUT_DIR}/bigwig/${sample}_FE.bw"
    
    bedGraphToBigWig "${OUTPUT_DIR}/peaks/${sample}_ppois.sorted.bdg" \
        "${CHROM_SIZES}" \
        "${OUTPUT_DIR}/bigwig/${sample}_ppois.bw"
    
    log "Fold enrichment processing complete for ${sample}."
}

compare_peaks() {
    local sample1=$1
    local sample2=$2
    local prefix=$3
    
    log "Comparing peaks between ${sample1} and ${sample2}..."
    
    # Create output directory
    mkdir -p "${OUTPUT_DIR}/differential/${prefix}"
    
    # Compare fold enrichment
    log "Comparing fold enrichment..."
    macs3 bdgdiff \
        --t1 "${OUTPUT_DIR}/peaks/${sample1}_FE.bdg" \
        --c1 "${OUTPUT_DIR}/peaks/${sample2}_FE.bdg" \
        --t2 "${OUTPUT_DIR}/peaks/${sample1}_FE.bdg" \
        --c2 "${OUTPUT_DIR}/peaks/${sample2}_FE.bdg" \
        --o-prefix "${prefix}" \
        --outdir "${OUTPUT_DIR}/differential/${prefix}_FE"
    
    # Compare p-values
    log "Comparing p-values..."
    macs3 bdgdiff \
        --t1 "${OUTPUT_DIR}/peaks/${sample1}_ppois.bdg" \
        --c1 "${OUTPUT_DIR}/peaks/${sample2}_ppois.bdg" \
        --t2 "${OUTPUT_DIR}/peaks/${sample1}_ppois.bdg" \
        --c2 "${OUTPUT_DIR}/peaks/${sample2}_ppois.bdg" \
        --o-prefix "${prefix}_ppois" \
        --outdir "${OUTPUT_DIR}/differential/${prefix}_ppois"
    
    log "Peak comparison complete for ${prefix}."
}

count_reads_in_peaks() {
    local sample=$1
    
    log "Creating SAF file for ${sample}..."
    awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' \
        "${OUTPUT_DIR}/peaks/${sample}/${sample}_peaks.broadPeak" > \
        "${OUTPUT_DIR}/peaks/${sample}_peaks.saf"
    
    log "Counting reads in peaks for ${sample}..."
    featureCounts -p -a "${OUTPUT_DIR}/peaks/${sample}_peaks.saf" -F SAF \
        -o "${OUTPUT_DIR}/peaks/${sample}-readCountInPeaks.txt" \
        "${OUTPUT_DIR}/filtered/${sample}.shifted.bam"
    
    log "Read counting complete for ${sample}."
}

# MAIN SCRIPT
##############################################################################

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--samples)
            SAMPLE_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            GENOME_INDEX="${HOME}/Databases/${GENOME}/${GENOME}_bt2index"
            BLACKLIST="${HOME}/Databases/${GENOME}/${GENOME}-blacklist.v2.bed"
            CHROM_SIZES="${HOME}/Databases/${GENOME}/${GENOME}.chrom.sizes"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
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

# Check if sample file exists
if [ ! -f "${SAMPLE_FILE}" ]; then
    log "ERROR: Sample file ${SAMPLE_FILE} not found."
    exit 1
fi

# Check dependencies
check_dependencies

# Create output directories
create_directories

# Process each sample
while read sample; do
    log "Processing sample: ${sample}"
    
    # Quality control
    run_fastqc "${sample}"
    
    # Trimming
    run_trimming "${sample}"
    
    # Alignment
    run_alignment "${sample}"
    
    # Filtering
    filter_mitochondrial "${sample}"
    remove_duplicates "${sample}"
    filter_blacklist "${sample}"
    shift_reads "${sample}"
    
    # Peak calling
    call_peaks "${sample}"
    
    # Generate bigWigs
    generate_bigwigs "${sample}"
    generate_enrichment "${sample}"
    
    # Count reads in peaks
    count_reads_in_peaks "${sample}"
    
    log "Processing complete for ${sample}"
done < "${SAMPLE_FILE}"

# Compare peaks between samples (example with Mut vs Ctrl)
# Uncomment and modify as needed for your specific comparisons
# compare_peaks "Mut-1" "Ctrl-1" "Mut_vs_Ctrl"

log "ATAC-seq pipeline completed successfully."
exit 0
