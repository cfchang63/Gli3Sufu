#!/usr/bin/env bash
#
# Split-DamID Analysis Workflow
# 
# This script automates the process of analyzing Split-DamID data 
# from raw fastq files to peak calling and coverage analysis.
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
GENOME_FASTA="${HOME}/Databases/${GENOME}/${GENOME}.fa"
BOWTIE2_INDEX="${HOME}/Databases/${GENOME}/${GENOME}_bt2index"
GATC_BED="${HOME}/Databases/${GENOME}/${GENOME}_GATC_sites.bed"
BLACKLIST="${HOME}/Databases/${GENOME}/blacklist_mtDNA.bed"
OUTPUT_DIR="./splitdamid_results"
SAMPLE_FILE="samples.txt"
BIN_SIZE=75
THREADS=24
GENOME_SIZE="mm" # mm for mouse, hs for human
EXTENSION_SCRIPT="${HOME}/Documents/Code/Python/SplitDamIDAnalysis/extend_reads_to_GATC_revised.py"
CALCULATE_LOG2FC=true
DAM_ONLY_PATTERN="DAM-[0-9]"
TREAT_PATTERN="Gli3-Hand2-[0-9]"
DOX_PATTERN="Dox"
MERGE_PEAKS=true
GENERATE_COUNTS=true

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -s, --samples FILE       Text file with sample names (default: ${SAMPLE_FILE})"
    echo "  -o, --output DIR         Output directory (default: ${OUTPUT_DIR})"
    echo "  -g, --genome NAME        Genome name (default: ${GENOME})"
    echo "  --fasta FILE             Genome FASTA file (default: ${GENOME_FASTA})"
    echo "  --bowtie2 DIR            Bowtie2 index directory (default: ${BOWTIE2_INDEX})"
    echo "  --gatc FILE              GATC sites BED file (default: ${GATC_BED})"
    echo "  --blacklist FILE         Blacklist+mtDNA regions BED file (default: ${BLACKLIST})"
    echo "  --extension-script FILE  Path to GATC extension script (default: ${EXTENSION_SCRIPT})"
    echo "  --bin-size NUM           Size of genome bins in bp (default: ${BIN_SIZE})"
    echo "  -t, --threads NUM        Number of threads (default: ${THREADS})"
    echo "  --genome-size SIZE       Genome size for MACS3 (mm or hs, default: ${GENOME_SIZE})"
    echo "  --control-pattern STR    Regex pattern to identify DAM-only control samples (default: ${DAM_ONLY_PATTERN})"
    echo "  --treat-pattern STR      Regex pattern to identify DAM-fusion treatment samples (default: ${TREAT_PATTERN})"
    echo "  --dox-pattern STR        Regex pattern to identify doxycycline-treated samples (default: ${DOX_PATTERN})"
    echo "  --skip-log2fc            Skip calculation of log2 fold change"
    echo "  --skip-merge-peaks       Skip merging of peaks across samples"
    echo "  --skip-counts            Skip generation of count matrix for DiffBind"
    echo "  -h, --help               Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    local dependencies=(bowtie2 samtools bedtools bamCoverage python3 conda)
    local missing=()
    
    log "Checking dependencies..."
    for cmd in "${dependencies[@]}"; do
        if ! command -v "${cmd}" &> /dev/null; then
            missing+=("${cmd}")
        fi
    done
    
    # Check if MACS3 is available in conda
    if ! conda env list | grep -q "MACS3"; then
        log "WARNING: MACS3 conda environment not found. Will create it."
        conda create -y -n MACS3 python=3.9 macs3
    fi
    
    if [ ${#missing[@]} -gt 0 ]; then
        log "ERROR: Missing dependencies: ${missing[*]}"
        log "Please install missing dependencies and try again."
        exit 1
    fi
    
    # Check for the extension script
    if [ ! -f "${EXTENSION_SCRIPT}" ]; then
        log "ERROR: GATC extension script not found at ${EXTENSION_SCRIPT}"
        exit 1
    fi
    
    log "All dependencies found."
}

create_directories() {
    local directories=(
        "${OUTPUT_DIR}"
        "${OUTPUT_DIR}/bam"
        "${OUTPUT_DIR}/bed"
        "${OUTPUT_DIR}/bigwig"
        "${OUTPUT_DIR}/peaks"
        "${OUTPUT_DIR}/bedgraph"
        "${OUTPUT_DIR}/analysis"
        "${OUTPUT_DIR}/log2fc"
        "${OUTPUT_DIR}/merged_peaks"
        "${OUTPUT_DIR}/counts"
    )
    
    log "Creating output directories..."
    for dir in "${directories[@]}"; do
        mkdir -p "${dir}"
    done
}

create_genome_bins() {
    log "Creating genome bins of size ${BIN_SIZE}bp..."
    
    # Check if genome FAI file exists
    if [ ! -f "${GENOME_FASTA}.fai" ]; then
        log "Generating FAI index for genome..."
        samtools faidx "${GENOME_FASTA}"
    fi
    
    # Create genome bins
    bedtools makewindows -g "${GENOME_FASTA}.fai" -w "${BIN_SIZE}" > "${OUTPUT_DIR}/genome_bins_${BIN_SIZE}bp.bed"
    
    log "Genome bins created."
}

process_sample() {
    local sample="$1"
    
    log "Processing sample: ${sample}"
    
    # Check if input files exist
    if [ ! -f "${sample}_R1.fastq.gz" ] || [ ! -f "${sample}_R2.fastq.gz" ]; then
        log "ERROR: Input files for ${sample} not found. Looking for ${sample}_R1.fastq.gz and ${sample}_R2.fastq.gz"
        return 1
    fi
    
    # Trimming with Trimmomatic (uncomment if needed)
    # log "Trimming adapters for ${sample}..."
    # java -jar trimmomatic.jar PE -threads ${THREADS} \
    #     "${sample}_R1.fastq.gz" "${sample}_R2.fastq.gz" \
    #     "${OUTPUT_DIR}/trimmed_${sample}_R1.fastq.gz" "${OUTPUT_DIR}/unpaired_${sample}_R1.fastq.gz" \
    #     "${OUTPUT_DIR}/trimmed_${sample}_R2.fastq.gz" "${OUTPUT_DIR}/unpaired_${sample}_R2.fastq.gz" \
    #     LEADING:5 TRAILING:5 MINLEN:50
    
    # For now we assume trimming is already done
    log "Assuming trimmed files exist as: trimmed_${sample}_R1.fastq.gz and trimmed_${sample}_R2.fastq.gz"
    
    # Alignment with Bowtie2
    log "Aligning ${sample} to ${GENOME}..."
    bowtie2 --local --very-sensitive --no-discordant --no-unal --threads "${THREADS}" -X 2000 \
        -x "${BOWTIE2_INDEX}" \
        -1 "trimmed_${sample}_R1.fastq.gz" -2 "trimmed_${sample}_R2.fastq.gz" \
        -S "${OUTPUT_DIR}/bam/${sample}_aligned.sam"
    
    # Convert SAM to BAM and sort by read name
    log "Converting SAM to BAM and sorting by read name for ${sample}..."
    samtools view -bS "${OUTPUT_DIR}/bam/${sample}_aligned.sam" | \
        samtools sort -n -o "${OUTPUT_DIR}/bam/${sample}_sorted.bam"
    
    # Remove SAM file to save space
    rm "${OUTPUT_DIR}/bam/${sample}_aligned.sam"
    
    # Convert BAM to BEDPE format
    log "Converting BAM to BEDPE for ${sample}..."
    bedtools bamtobed -bedpe -i "${OUTPUT_DIR}/bam/${sample}_sorted.bam" > "${OUTPUT_DIR}/bed/${sample}_pe.bed"
    
    # Extend reads to the next GATC site
    log "Extending reads to GATC sites for ${sample}..."
    python3 "${EXTENSION_SCRIPT}" \
        "${OUTPUT_DIR}/bed/${sample}_pe.bed" "${GATC_BED}" "${OUTPUT_DIR}/bed/${sample}_extended.bed"
    
    # Convert extended BEDPE back to BAM
    log "Converting extended BEDPE back to BAM for ${sample}..."
    bedtools bedpetobam -i "${OUTPUT_DIR}/bed/${sample}_extended.bed" -g "${GENOME_FASTA}.fai" > \
        "${OUTPUT_DIR}/bam/${sample}_extended.bam"
    
    # Sort and index extended BAM
    log "Sorting and indexing extended BAM for ${sample}..."
    samtools sort -o "${OUTPUT_DIR}/bam/${sample}_extended_sorted.bam" "${OUTPUT_DIR}/bam/${sample}_extended.bam"
    samtools index "${OUTPUT_DIR}/bam/${sample}_extended_sorted.bam"
    
    # Remove intermediate file
    rm "${OUTPUT_DIR}/bam/${sample}_extended.bam"
    
    # Filter blacklisted regions and mtDNA
    log "Filtering blacklisted regions and mtDNA for ${sample}..."
    bedtools intersect -v -abam "${OUTPUT_DIR}/bam/${sample}_extended_sorted.bam" -b "${BLACKLIST}" > \
        "${OUTPUT_DIR}/bam/${sample}_filtered.bam"
    
    # Sort and index filtered BAM
    log "Sorting and indexing filtered BAM for ${sample}..."
    samtools sort -o "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" "${OUTPUT_DIR}/bam/${sample}_filtered.bam"
    samtools index "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam"
    
    # Remove intermediate file
    rm "${OUTPUT_DIR}/bam/${sample}_filtered.bam"
    
    # Generate BigWig file
    log "Generating BigWig for ${sample}..."
    bamCoverage -b "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" \
        -o "${OUTPUT_DIR}/bigwig/${sample}.bw" \
        --normalizeUsing CPM \
        --numberOfProcessors "${THREADS}" \
        --blackListFileName "${BLACKLIST}"
    
    # Generate non-binned genome coverage
    log "Generating non-binned coverage for ${sample}..."
    bedtools genomecov -ibam "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" -bg > \
        "${OUTPUT_DIR}/bedgraph/${sample}_nonbinned.bedgraph"
    
    # Generate binned genome coverage
    log "Counting reads in bins for ${sample}..."
    bedtools coverage -a "${OUTPUT_DIR}/genome_bins_${BIN_SIZE}bp.bed" \
        -b "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" -counts > \
        "${OUTPUT_DIR}/analysis/${sample}_bin_counts.bed"
    
    # Call peaks with MACS3
    log "Calling peaks for ${sample}..."
    # Set genome size parameter
    if [ "${GENOME_SIZE}" = "mm" ]; then
        local macs_genome="2652783500"
    elif [ "${GENOME_SIZE}" = "hs" ]; then
        local macs_genome="3209286100"
    else
        local macs_genome="${GENOME_SIZE}"
    fi
    
    # Activate MACS3 conda environment
    conda activate MACS3
    
    # Call peaks
    macs3 callpeak -t "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" \
        -f BAMPE -g "${macs_genome}" --broad -B \
        -n "${sample}_broad" --outdir "${OUTPUT_DIR}/peaks/${sample}" \
        --keep-dup all
    
    # Deactivate MACS3 environment
    conda deactivate
    
    log "Processing completed for ${sample}"
}

merge_peaks_across_samples() {
    log "Merging peaks across samples..."
    
    # Get all sample names
    local all_samples=$(cat "${SAMPLE_FILE}")
    
    # Create arrays for different groups
    local dam_samples=()
    local dam_dox_samples=()
    local treat_samples=()
    local treat_dox_samples=()
    
    # Categorize samples based on patterns
    for sample in ${all_samples}; do
        if [[ "${sample}" =~ ${DAM_ONLY_PATTERN} && ! "${sample}" =~ ${DOX_PATTERN} ]]; then
            dam_samples+=("${sample}")
        elif [[ "${sample}" =~ ${DAM_ONLY_PATTERN} && "${sample}" =~ ${DOX_PATTERN} ]]; then
            dam_dox_samples+=("${sample}")
        elif [[ "${sample}" =~ ${TREAT_PATTERN} && ! "${sample}" =~ ${DOX_PATTERN} ]]; then
            treat_samples+=("${sample}")
        elif [[ "${sample}" =~ ${TREAT_PATTERN} && "${sample}" =~ ${DOX_PATTERN} ]]; then
            treat_dox_samples+=("${sample}")
        fi
    done
    
    # Define groups
    local groups=("DAM" "DAM-Dox" "Gli3-Hand2" "Gli3-Hand2-Dox")
    
    # Process each group
    for group in "${groups[@]}"; do
        log "Merging peaks for group: ${group}"
        
        # Get list of peak files for the group
        local peak_files=""
        
        if [ "${group}" = "DAM" ]; then
            for sample in "${dam_samples[@]}"; do
                peak_files="${peak_files} ${OUTPUT_DIR}/peaks/${sample}/${sample}_broad_peaks.broadPeak"
            done
        elif [ "${group}" = "DAM-Dox" ]; then
            for sample in "${dam_dox_samples[@]}"; do
                peak_files="${peak_files} ${OUTPUT_DIR}/peaks/${sample}/${sample}_broad_peaks.broadPeak"
            done
        elif [ "${group}" = "Gli3-Hand2" ]; then
            for sample in "${treat_samples[@]}"; do
                peak_files="${peak_files} ${OUTPUT_DIR}/peaks/${sample}/${sample}_broad_peaks.broadPeak"
            done
        elif [ "${group}" = "Gli3-Hand2-Dox" ]; then
            for sample in "${treat_dox_samples[@]}"; do
                peak_files="${peak_files} ${OUTPUT_DIR}/peaks/${sample}/${sample}_broad_peaks.broadPeak"
            done
        fi
        
        # Skip if no peak files found
        if [ -z "${peak_files}" ]; then
            log "WARNING: No peak files found for group ${group}"
            continue
        fi
        
        # Concatenate peak files
        cat ${peak_files} > "${OUTPUT_DIR}/merged_peaks/${group}_combined_peaks.bed"
        
        # Sort and merge overlapping peaks
        sort -k1,1 -k2,2n "${OUTPUT_DIR}/merged_peaks/${group}_combined_peaks.bed" | \
            bedtools merge > "${OUTPUT_DIR}/merged_peaks/${group}_merged_peaks.bed"
        
        log "Created merged peaks for ${group}: ${OUTPUT_DIR}/merged_peaks/${group}_merged_peaks.bed"
    done
    
    # Create a unified peak set by merging peaks across all groups
    log "Creating unified peak set across all groups..."
    
    # Check if all group peak files exist
    for group in "${groups[@]}"; do
        if [ ! -f "${OUTPUT_DIR}/merged_peaks/${group}_merged_peaks.bed" ]; then
            log "WARNING: Missing merged peaks for group ${group}"
        fi
    done
    
    # Concatenate all group peak files
    cat "${OUTPUT_DIR}/merged_peaks/DAM_merged_peaks.bed" \
        "${OUTPUT_DIR}/merged_peaks/DAM-Dox_merged_peaks.bed" \
        "${OUTPUT_DIR}/merged_peaks/Gli3-Hand2_merged_peaks.bed" \
        "${OUTPUT_DIR}/merged_peaks/Gli3-Hand2-Dox_merged_peaks.bed" > \
        "${OUTPUT_DIR}/merged_peaks/all_peaks.bed"
    
    # Sort and merge to create unified peak set
    sort -k1,1 -k2,2n "${OUTPUT_DIR}/merged_peaks/all_peaks.bed" | \
        bedtools merge > "${OUTPUT_DIR}/merged_peaks/unified_peaks.bed"
    
    log "Created unified peak set: ${OUTPUT_DIR}/merged_peaks/unified_peaks.bed"
}

generate_count_matrix() {
    log "Generating count matrix for differential analysis..."
    
    # Check if unified peak set exists
    if [ ! -f "${OUTPUT_DIR}/merged_peaks/unified_peaks.bed" ]; then
        log "ERROR: Unified peak set not found. Run merge_peaks_across_samples first."
        return 1
    fi
    
    # Get all sample names
    local all_samples=$(cat "${SAMPLE_FILE}")
    
    # Count for each sample
    for sample in ${all_samples}; do
        log "Counting reads for sample: ${sample}"
        
        # Check if filtered BAM exists
        if [ ! -f "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" ]; then
            log "WARNING: Filtered BAM file not found for ${sample}, using filtered.bam if available"
            if [ -f "${OUTPUT_DIR}/bam/${sample}_filtered.bam" ]; then
                bedtools coverage -a "${OUTPUT_DIR}/merged_peaks/unified_peaks.bed" \
                    -b "${OUTPUT_DIR}/bam/${sample}_filtered.bam" -counts > \
                    "${OUTPUT_DIR}/counts/${sample}_counts.txt"
            else
                log "WARNING: No suitable BAM file found for ${sample}, skipping"
                continue
            fi
        else
            bedtools coverage -a "${OUTPUT_DIR}/merged_peaks/unified_peaks.bed" \
                -b "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" -counts > \
                "${OUTPUT_DIR}/counts/${sample}_counts.txt"
        fi
    done
    
    # Create combined count matrix using Python
    log "Creating combined count matrix..."
    
    # Create Python script for combining counts
    cat > "${OUTPUT_DIR}/counts/combine_counts.py" << 'EOF'
#!/usr/bin/env python3
import os
import glob
import pandas as pd
import sys

# Directory containing count files
counts_dir = sys.argv[1]
output_file = sys.argv[2]

# Get all count files
count_files = glob.glob(os.path.join(counts_dir, "*_counts.txt"))

if not count_files:
    print("No count files found!")
    sys.exit(1)

# Initialize DataFrame with coordinates from first file
df = pd.read_csv(count_files[0], sep='\t', header=None, 
                names=['Chr', 'Start', 'End', count_files[0].split('/')[-1].replace('_counts.txt', '')])

# Set index for joining
df.set_index(['Chr', 'Start', 'End'], inplace=True)

# Add counts from other files
for file in count_files[1:]:
    sample_name = file.split('/')[-1].replace('_counts.txt', '')
    temp_df = pd.read_csv(file, sep='\t', header=None, 
                        names=['Chr', 'Start', 'End', sample_name])
    temp_df.set_index(['Chr', 'Start', 'End'], inplace=True)
    df = df.join(temp_df)

# Reset index to get Chr, Start, End as columns
df.reset_index(inplace=True)

# Save combined matrix
df.to_csv(output_file, sep='\t', index=False)
print(f"Created combined count matrix with {len(df)} rows and {len(df.columns)} columns")
EOF
    
    # Make Python script executable
    chmod +x "${OUTPUT_DIR}/counts/combine_counts.py"
    
    # Run Python script
    python3 "${OUTPUT_DIR}/counts/combine_counts.py" \
        "${OUTPUT_DIR}/counts" \
        "${OUTPUT_DIR}/counts/combined_counts.tsv"
    
    log "Count matrix generation complete: ${OUTPUT_DIR}/counts/combined_counts.tsv"
}

calculate_log2fc() {
    log "Calculating log2 fold changes..."
    
    # Create directory for pairwise comparisons
    mkdir -p "${OUTPUT_DIR}/log2fc"
    
    # Get all sample names
    local all_samples=$(cat "${SAMPLE_FILE}")
    
    # Create arrays for control and treatment samples
    local dam_samples=()
    local dam_dox_samples=()
    local treat_samples=()
    local treat_dox_samples=()
    
    # Categorize samples based on patterns
    for sample in ${all_samples}; do
        if [[ "${sample}" =~ ${DAM_ONLY_PATTERN} && ! "${sample}" =~ ${DOX_PATTERN} ]]; then
            dam_samples+=("${sample}")
        elif [[ "${sample}" =~ ${DAM_ONLY_PATTERN} && "${sample}" =~ ${DOX_PATTERN} ]]; then
            dam_dox_samples+=("${sample}")
        elif [[ "${sample}" =~ ${TREAT_PATTERN} && ! "${sample}" =~ ${DOX_PATTERN} ]]; then
            treat_samples+=("${sample}")
        elif [[ "${sample}" =~ ${TREAT_PATTERN} && "${sample}" =~ ${DOX_PATTERN} ]]; then
            treat_dox_samples+=("${sample}")
        fi
    done
    
    log "DAM controls: ${dam_samples[*]}"
    log "DAM+DOX controls: ${dam_dox_samples[*]}"
    log "Treatment samples: ${treat_samples[*]}"
    log "Treatment+DOX samples: ${treat_dox_samples[*]}"
    
    # Calculate log2fc for each treatment vs corresponding control
    for treat in "${treat_samples[@]}"; do
        for dam in "${dam_samples[@]}"; do
            log "Calculating log2FC for ${treat} vs ${dam}..."
            
            # Create R script for log2fc calculation
            cat > "${OUTPUT_DIR}/log2fc/calculate_log2fc_${treat}_vs_${dam}.R" << EOF
# Log2FC calculation for ${treat} vs ${dam}
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)

# Load binned data
treat_bins <- read.table("${OUTPUT_DIR}/analysis/${treat}_bin_counts.bed", 
                         col.names=c("chr", "start", "end", "count"))
control_bins <- read.table("${OUTPUT_DIR}/analysis/${dam}_bin_counts.bed", 
                           col.names=c("chr", "start", "end", "count"))

# Add pseudocount and calculate log2FC
treat_bins\$count_norm <- treat_bins\$count + 1
control_bins\$count_norm <- control_bins\$count + 1

# Merge data
merged_data <- data.frame(
  chr = treat_bins\$chr,
  start = treat_bins\$start,
  end = treat_bins\$end,
  treat_count = treat_bins\$count,
  control_count = control_bins\$count,
  log2FC = log2(treat_bins\$count_norm / control_bins\$count_norm)
)

# Write output
write.table(merged_data, file="${OUTPUT_DIR}/log2fc/${treat}_vs_${dam}_log2fc.txt",
            quote=FALSE, sep="\t", row.names=FALSE)

# Create bedGraph for visualization
bedgraph_data <- data.frame(
  chr = merged_data\$chr,
  start = merged_data\$start,
  end = merged_data\$end,
  score = merged_data\$log2FC
)

# Filter out extreme values
bedgraph_data\$score[bedgraph_data\$score > 5] <- 5
bedgraph_data\$score[bedgraph_data\$score < -5] <- -5

# Write bedGraph
write.table(bedgraph_data, file="${OUTPUT_DIR}/log2fc/${treat}_vs_${dam}_log2fc.bedgraph",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Create BigWig if possible
if (require(rtracklayer)) {
  try({
    gr <- GRanges(
      seqnames = bedgraph_data\$chr,
      ranges = IRanges(start = bedgraph_data\$start + 1, end = bedgraph_data\$end),
      score = bedgraph_data\$score
    )
    export(gr, "${OUTPUT_DIR}/log2fc/${treat}_vs_${dam}_log2fc.bw")
  })
}

cat("Log2FC calculation complete for ${treat} vs ${dam}\\n")
EOF
            
            # Run R script
            Rscript "${OUTPUT_DIR}/log2fc/calculate_log2fc_${treat}_vs_${dam}.R"
        done
    done
    
    # Similarly for DOX-treated samples
    for treat_dox in "${treat_dox_samples[@]}"; do
        for dam_dox in "${dam_dox_samples[@]}"; do
            log "Calculating log2FC for ${treat_dox} vs ${dam_dox}..."
            
            # Create R script for log2fc calculation
            cat > "${OUTPUT_DIR}/log2fc/calculate_log2fc_${treat_dox}_vs_${dam_dox}.R" << EOF
# Log2FC calculation for ${treat_dox} vs ${dam_dox}
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)

# Load binned data
treat_bins <- read.table("${OUTPUT_DIR}/analysis/${treat_dox}_bin_counts.bed", 
                         col.names=c("chr", "start", "end", "count"))
control_bins <- read.table("${OUTPUT_DIR}/analysis/${dam_dox}_bin_counts.bed", 
                           col.names=c("chr", "start", "end", "count"))

# Add pseudocount and calculate log2FC
treat_bins\$count_norm <- treat_bins\$count + 1
control_bins\$count_norm <- control_bins\$count + 1

# Merge data
merged_data <- data.frame(
  chr = treat_bins\$chr,
  start = treat_bins\$start,
  end = treat_bins\$end,
  treat_count = treat_bins\$count,
  control_count = control_bins\$count,
  log2FC = log2(treat_bins\$count_norm / control_bins\$count_norm)
)

# Write output
write.table(merged_data, file="${OUTPUT_DIR}/log2fc/${treat_dox}_vs_${dam_dox}_log2fc.txt",
            quote=FALSE, sep="\t", row.names=FALSE)

# Create bedGraph for visualization
bedgraph_data <- data.frame(
  chr = merged_data\$chr,
  start = merged_data\$start,
  end = merged_data\$end,
  score = merged_data\$log2FC
)

# Filter out extreme values
bedgraph_data\$score[bedgraph_data\$score > 5] <- 5
bedgraph_data\$score[bedgraph_data\$score < -5] <- -5

# Write bedGraph
write.table(bedgraph_data, file="${OUTPUT_DIR}/log2fc/${treat_dox}_vs_${dam_dox}_log2fc.bedgraph",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Create BigWig if possible
if (require(rtracklayer)) {
  try({
    gr <- GRanges(
      seqnames = bedgraph_data\$chr,
      ranges = IRanges(start = bedgraph_data\$start + 1, end = bedgraph_data\$end),
      score = bedgraph_data\$score
    )
    export(gr, "${OUTPUT_DIR}/log2fc/${treat_dox}_vs_${dam_dox}_log2fc.bw")
  })
}

cat("Log2FC calculation complete for ${treat_dox} vs ${dam_dox}\\n")
EOF
            
            # Run R script
            Rscript "${OUTPUT_DIR}/log2fc/calculate_log2fc_${treat_dox}_vs_${dam_dox}.R"
        done
    done
    
    log "Log2 fold change calculations complete."
}

create_sample_file() {
    log "No sample file provided. Creating template at ${SAMPLE_FILE}..."
    
    cat > "${SAMPLE_FILE}" << EOF
DAM-1
DAM-2
DAM-3
DAM-Dox-1
DAM-Dox-2
DAM-Dox-3
Gli3-Hand2-1
Gli3-Hand2-2
Gli3-Hand2-3
Gli3-Hand2-Dox-1
Gli3-Hand2-Dox-2
Gli3-Hand2-Dox-3
EOF
    
    log "Sample file template created. Edit this file with your sample names."
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
            shift 2
            ;;
        --fasta)
            GENOME_FASTA="$2"
            shift 2
            ;;
        --bowtie2)
            BOWTIE2_INDEX="$2"
            shift 2
            ;;
        --gatc)
            GATC_BED="$2"
            shift 2
            ;;
        --blacklist)
            BLACKLIST="$2"
            shift 2
            ;;
        --extension-script)
            EXTENSION_SCRIPT="$2"
            shift 2
            ;;
        --bin-size)
            BIN_SIZE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --genome-size)
            GENOME_SIZE="$2"
            shift 2
            ;;
        --control-pattern)
            DAM_ONLY_PATTERN="$2"
            shift 2
            ;;
        --treat-pattern)
            TREAT_PATTERN="$2"
            shift 2
            ;;
        --dox-pattern)
            DOX_PATTERN="$2"
            shift 2
            ;;
        --skip-log2fc)
            CALCULATE_LOG2FC=false
            shift
            ;;
        --skip-merge-peaks)
            MERGE_PEAKS=false
            shift
            ;;
        --skip-counts)
            GENERATE_COUNTS=false
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

# Check if sample file exists, create template if not
if [ ! -f "${SAMPLE_FILE}" ]; then
    create_sample_file
    log "Please edit the sample file and run the script again."
    exit 1
fi

# Check dependencies
check_dependencies

# Create output directories
create_directories

# Create genome bins
create_genome_bins

# Process each sample
while read sample; do
    process_sample "${sample}"
done < "${SAMPLE_FILE}"

# Merge peaks across samples if requested
if [ "${MERGE_PEAKS}" = true ]; then
    merge_peaks_across_samples
    
    # Generate count matrix if requested
    if [ "${GENERATE_COUNTS}" = true ]; then
        generate_count_matrix
    fi
fi

# Calculate log2 fold changes if requested
if [ "${CALCULATE_LOG2FC}" = true ]; then
    calculate_log2fc
fi

log "Split-DamID analysis workflow completed successfully!"
log "Results can be found in: ${OUTPUT_DIR}"

# Print summary
echo ""
echo "Summary of results:"
echo "==================="
echo "Processed $(wc -l < ${SAMPLE_FILE}) samples"
echo "BAM files: ${OUTPUT_DIR}/bam/"
echo "BigWig files: ${OUTPUT_DIR}/bigwig/"
echo "Binned count data: ${OUTPUT_DIR}/analysis/"
echo "Peak calls: ${OUTPUT_DIR}/peaks/"

if [ "${MERGE_PEAKS}" = true ]; then
    echo "Merged peaks: ${OUTPUT_DIR}/merged_peaks/"
    if [ "${GENERATE_COUNTS}" = true ]; then
        echo "Count matrix: ${OUTPUT_DIR}/counts/combined_counts.tsv"
    fi
fi

if [ "${CALCULATE_LOG2FC}" = true ]; then
    echo "Log2 fold change data: ${OUTPUT_DIR}/log2fc/"
fi

echo ""
log "Done."
exit 0#!/usr/bin/env bash
#
# Split-DamID Analysis Workflow
# 
# This script automates the process of analyzing Split-DamID data 
# from raw fastq files to peak calling and coverage analysis.
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
GENOME_FASTA="${HOME}/Databases/${GENOME}/${GENOME}.fa"
BOWTIE2_INDEX="${HOME}/Databases/${GENOME}/${GENOME}_bt2index"
GATC_BED="${HOME}/Databases/${GENOME}/${GENOME}_GATC_sites.bed"
BLACKLIST="${HOME}/Databases/${GENOME}/blacklist_mtDNA.bed"
OUTPUT_DIR="./splitdamid_results"
SAMPLE_FILE="samples.txt"
BIN_SIZE=75
THREADS=24
GENOME_SIZE="mm" # mm for mouse, hs for human
EXTENSION_SCRIPT="${HOME}/Documents/Code/Python/SplitDamIDAnalysis/extend_reads_to_GATC_revised.py"
CALCULATE_LOG2FC=true
DAM_ONLY_PATTERN="DAM-[0-9]"
TREAT_PATTERN="Gli3-Hand2-[0-9]"
DOX_PATTERN="Dox"

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -s, --samples FILE       Text file with sample names (default: ${SAMPLE_FILE})"
    echo "  -o, --output DIR         Output directory (default: ${OUTPUT_DIR})"
    echo "  -g, --genome NAME        Genome name (default: ${GENOME})"
    echo "  --fasta FILE             Genome FASTA file (default: ${GENOME_FASTA})"
    echo "  --bowtie2 DIR            Bowtie2 index directory (default: ${BOWTIE2_INDEX})"
    echo "  --gatc FILE              GATC sites BED file (default: ${GATC_BED})"
    echo "  --blacklist FILE         Blacklist+mtDNA regions BED file (default: ${BLACKLIST})"
    echo "  --extension-script FILE  Path to GATC extension script (default: ${EXTENSION_SCRIPT})"
    echo "  --bin-size NUM           Size of genome bins in bp (default: ${BIN_SIZE})"
    echo "  -t, --threads NUM        Number of threads (default: ${THREADS})"
    echo "  --genome-size SIZE       Genome size for MACS3 (mm or hs, default: ${GENOME_SIZE})"
    echo "  --control-pattern STR    Regex pattern to identify DAM-only control samples (default: ${DAM_ONLY_PATTERN})"
    echo "  --treat-pattern STR      Regex pattern to identify DAM-fusion treatment samples (default: ${TREAT_PATTERN})"
    echo "  --dox-pattern STR        Regex pattern to identify doxycycline-treated samples (default: ${DOX_PATTERN})"
    echo "  --skip-log2fc            Skip calculation of log2 fold change"
    echo "  -h, --help               Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    local dependencies=(bowtie2 samtools bedtools bamCoverage python3 conda)
    local missing=()
    
    log "Checking dependencies..."
    for cmd in "${dependencies[@]}"; do
        if ! command -v "${cmd}" &> /dev/null; then
            missing+=("${cmd}")
        fi
    done
    
    # Check if MACS3 is available in conda
    if ! conda env list | grep -q "MACS3"; then
        log "WARNING: MACS3 conda environment not found. Will create it."
        conda create -y -n MACS3 python=3.9 macs3
    fi
    
    if [ ${#missing[@]} -gt 0 ]; then
        log "ERROR: Missing dependencies: ${missing[*]}"
        log "Please install missing dependencies and try again."
        exit 1
    fi
    
    # Check for the extension script
    if [ ! -f "${EXTENSION_SCRIPT}" ]; then
        log "ERROR: GATC extension script not found at ${EXTENSION_SCRIPT}"
        exit 1
    fi
    
    log "All dependencies found."
}

create_directories() {
    local directories=(
        "${OUTPUT_DIR}"
        "${OUTPUT_DIR}/bam"
        "${OUTPUT_DIR}/bed"
        "${OUTPUT_DIR}/bigwig"
        "${OUTPUT_DIR}/peaks"
        "${OUTPUT_DIR}/bedgraph"
        "${OUTPUT_DIR}/analysis"
        "${OUTPUT_DIR}/log2fc"
    )
    
    log "Creating output directories..."
    for dir in "${directories[@]}"; do
        mkdir -p "${dir}"
    done
}

create_genome_bins() {
    log "Creating genome bins of size ${BIN_SIZE}bp..."
    
    # Check if genome FAI file exists
    if [ ! -f "${GENOME_FASTA}.fai" ]; then
        log "Generating FAI index for genome..."
        samtools faidx "${GENOME_FASTA}"
    fi
    
    # Create genome bins
    bedtools makewindows -g "${GENOME_FASTA}.fai" -w "${BIN_SIZE}" > "${OUTPUT_DIR}/genome_bins_${BIN_SIZE}bp.bed"
    
    log "Genome bins created."
}

process_sample() {
    local sample="$1"
    
    log "Processing sample: ${sample}"
    
    # Check if input files exist
    if [ ! -f "${sample}_R1.fastq.gz" ] || [ ! -f "${sample}_R2.fastq.gz" ]; then
        log "ERROR: Input files for ${sample} not found. Looking for ${sample}_R1.fastq.gz and ${sample}_R2.fastq.gz"
        return 1
    fi
    
    # Trimming with Trimmomatic (uncomment if needed)
    # log "Trimming adapters for ${sample}..."
    # java -jar trimmomatic.jar PE -threads ${THREADS} \
    #     "${sample}_R1.fastq.gz" "${sample}_R2.fastq.gz" \
    #     "${OUTPUT_DIR}/trimmed_${sample}_R1.fastq.gz" "${OUTPUT_DIR}/unpaired_${sample}_R1.fastq.gz" \
    #     "${OUTPUT_DIR}/trimmed_${sample}_R2.fastq.gz" "${OUTPUT_DIR}/unpaired_${sample}_R2.fastq.gz" \
    #     LEADING:5 TRAILING:5 MINLEN:50
    
    # For now we assume trimming is already done
    log "Assuming trimmed files exist as: trimmed_${sample}_R1.fastq.gz and trimmed_${sample}_R2.fastq.gz"
    
    # Alignment with Bowtie2
    log "Aligning ${sample} to ${GENOME}..."
    bowtie2 --local --very-sensitive --no-discordant --no-unal --threads "${THREADS}" -X 2000 \
        -x "${BOWTIE2_INDEX}" \
        -1 "trimmed_${sample}_R1.fastq.gz" -2 "trimmed_${sample}_R2.fastq.gz" \
        -S "${OUTPUT_DIR}/bam/${sample}_aligned.sam"
    
    # Convert SAM to BAM and sort by read name
    log "Converting SAM to BAM and sorting by read name for ${sample}..."
    samtools view -bS "${OUTPUT_DIR}/bam/${sample}_aligned.sam" | \
        samtools sort -n -o "${OUTPUT_DIR}/bam/${sample}_sorted.bam"
    
    # Remove SAM file to save space
    rm "${OUTPUT_DIR}/bam/${sample}_aligned.sam"
    
    # Convert BAM to BEDPE format
    log "Converting BAM to BEDPE for ${sample}..."
    bedtools bamtobed -bedpe -i "${OUTPUT_DIR}/bam/${sample}_sorted.bam" > "${OUTPUT_DIR}/bed/${sample}_pe.bed"
    
    # Extend reads to the next GATC site
    log "Extending reads to GATC sites for ${sample}..."
    python3 "${EXTENSION_SCRIPT}" \
        "${OUTPUT_DIR}/bed/${sample}_pe.bed" "${GATC_BED}" "${OUTPUT_DIR}/bed/${sample}_extended.bed"
    
    # Convert extended BEDPE back to BAM
    log "Converting extended BEDPE back to BAM for ${sample}..."
    bedtools bedpetobam -i "${OUTPUT_DIR}/bed/${sample}_extended.bed" -g "${GENOME_FASTA}.fai" > \
        "${OUTPUT_DIR}/bam/${sample}_extended.bam"
    
    # Sort and index extended BAM
    log "Sorting and indexing extended BAM for ${sample}..."
    samtools sort -o "${OUTPUT_DIR}/bam/${sample}_extended_sorted.bam" "${OUTPUT_DIR}/bam/${sample}_extended.bam"
    samtools index "${OUTPUT_DIR}/bam/${sample}_extended_sorted.bam"
    
    # Remove intermediate file
    rm "${OUTPUT_DIR}/bam/${sample}_extended.bam"
    
    # Filter blacklisted regions and mtDNA
    log "Filtering blacklisted regions and mtDNA for ${sample}..."
    bedtools intersect -v -abam "${OUTPUT_DIR}/bam/${sample}_extended_sorted.bam" -b "${BLACKLIST}" > \
        "${OUTPUT_DIR}/bam/${sample}_filtered.bam"
    
    # Sort and index filtered BAM
    log "Sorting and indexing filtered BAM for ${sample}..."
    samtools sort -o "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" "${OUTPUT_DIR}/bam/${sample}_filtered.bam"
    samtools index "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam"
    
    # Remove intermediate file
    rm "${OUTPUT_DIR}/bam/${sample}_filtered.bam"
    
    # Generate BigWig file
    log "Generating BigWig for ${sample}..."
    bamCoverage -b "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" \
        -o "${OUTPUT_DIR}/bigwig/${sample}.bw" \
        --normalizeUsing CPM \
        --numberOfProcessors "${THREADS}" \
        --blackListFileName "${BLACKLIST}"
    
    # Generate non-binned genome coverage
    log "Generating non-binned coverage for ${sample}..."
    bedtools genomecov -ibam "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" -bg > \
        "${OUTPUT_DIR}/bedgraph/${sample}_nonbinned.bedgraph"
    
    # Generate binned genome coverage
    log "Counting reads in bins for ${sample}..."
    bedtools coverage -a "${OUTPUT_DIR}/genome_bins_${BIN_SIZE}bp.bed" \
        -b "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" -counts > \
        "${OUTPUT_DIR}/analysis/${sample}_bin_counts.bed"
    
    # Call peaks with MACS3
    log "Calling peaks for ${sample}..."
    # Set genome size parameter
    if [ "${GENOME_SIZE}" = "mm" ]; then
        local macs_genome="2652783500"
    elif [ "${GENOME_SIZE}" = "hs" ]; then
        local macs_genome="3209286100"
    else
        local macs_genome="${GENOME_SIZE}"
    fi
    
    # Activate MACS3 conda environment
    conda activate MACS3
    
    # Call peaks
    macs3 callpeak -t "${OUTPUT_DIR}/bam/${sample}_filtered_sorted.bam" \
        -f BAMPE -g "${macs_genome}" --broad -B \
        -n "${sample}_broad" --outdir "${OUTPUT_DIR}/peaks/${sample}" \
        --keep-dup all
    
    # Deactivate MACS3 environment
    conda deactivate
    
    log "Processing completed for ${sample}"
}

calculate_log2fc() {
    log "Calculating log2 fold changes..."
    
    # Create directory for pairwise comparisons
    mkdir -p "${OUTPUT_DIR}/log2fc"
    
    # Get all sample names
    local all_samples=$(cat "${SAMPLE_FILE}")
    
    # Create arrays for control and treatment samples
    local dam_samples=()
    local dam_dox_samples=()
    local treat_samples=()
    local treat_dox_samples=()
    
    # Categorize samples based on patterns
    for sample in ${all_samples}; do
        if [[ "${sample}" =~ ${DAM_ONLY_PATTERN} && ! "${sample}" =~ ${DOX_PATTERN} ]]; then
            dam_samples+=("${sample}")
        elif [[ "${sample}" =~ ${DAM_ONLY_PATTERN} && "${sample}" =~ ${DOX_PATTERN} ]]; then
            dam_dox_samples+=("${sample}")
        elif [[ "${sample}" =~ ${TREAT_PATTERN} && ! "${sample}" =~ ${DOX_PATTERN} ]]; then
            treat_samples+=("${sample}")
        elif [[ "${sample}" =~ ${TREAT_PATTERN} && "${sample}" =~ ${DOX_PATTERN} ]]; then
            treat_dox_samples+=("${sample}")
        fi
    done
    
    log "DAM controls: ${dam_samples[*]}"
    log "DAM+DOX controls: ${dam_dox_samples[*]}"
    log "Treatment samples: ${treat_samples[*]}"
    log "Treatment+DOX samples: ${treat_dox_samples[*]}"
    
    # Calculate log2fc for each treatment vs corresponding control
    for treat in "${treat_samples[@]}"; do
        for dam in "${dam_samples[@]}"; do
            log "Calculating log2FC for ${treat} vs ${dam}..."
            
            # Create R script for log2fc calculation
            cat > "${OUTPUT_DIR}/log2fc/calculate_log2fc_${treat}_vs_${dam}.R" << EOF
# Log2FC calculation for ${treat} vs ${dam}
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)

# Load binned data
treat_bins <- read.table("${OUTPUT_DIR}/analysis/${treat}_bin_counts.bed", 
                         col.names=c("chr", "start", "end", "count"))
control_bins <- read.table("${OUTPUT_DIR}/analysis/${dam}_bin_counts.bed", 
                           col.names=c("chr", "start", "end", "count"))

# Add pseudocount and calculate log2FC
treat_bins\$count_norm <- treat_bins\$count + 1
control_bins\$count_norm <- control_bins\$count + 1

# Merge data
merged_data <- data.frame(
  chr = treat_bins\$chr,
  start = treat_bins\$start,
  end = treat_bins\$end,
  treat_count = treat_bins\$count,
  control_count = control_bins\$count,
  log2FC = log2(treat_bins\$count_norm / control_bins\$count_norm)
)

# Write output
write.table(merged_data, file="${OUTPUT_DIR}/log2fc/${treat}_vs_${dam}_log2fc.txt",
            quote=FALSE, sep="\t", row.names=FALSE)

# Create bedGraph for visualization
bedgraph_data <- data.frame(
  chr = merged_data\$chr,
  start = merged_data\$start,
  end = merged_data\$end,
  score = merged_data\$log2FC
)

# Filter out extreme values
bedgraph_data\$score[bedgraph_data\$score > 5] <- 5
bedgraph_data\$score[bedgraph_data\$score < -5] <- -5

# Write bedGraph
write.table(bedgraph_data, file="${OUTPUT_DIR}/log2fc/${treat}_vs_${dam}_log2fc.bedgraph",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Create BigWig if possible
if (require(rtracklayer)) {
  try({
    gr <- GRanges(
      seqnames = bedgraph_data\$chr,
      ranges = IRanges(start = bedgraph_data\$start + 1, end = bedgraph_data\$end),
      score = bedgraph_data\$score
    )
    export(gr, "${OUTPUT_DIR}/log2fc/${treat}_vs_${dam}_log2fc.bw")
  })
}

cat("Log2FC calculation complete for ${treat} vs ${dam}\\n")
EOF
            
            # Run R script
            Rscript "${OUTPUT_DIR}/log2fc/calculate_log2fc_${treat}_vs_${dam}.R"
        done
    done
    
    # Similarly for DOX-treated samples
    for treat_dox in "${treat_dox_samples[@]}"; do
        for dam_dox in "${dam_dox_samples[@]}"; do
            log "Calculating log2FC for ${treat_dox} vs ${dam_dox}..."
            
            # Create R script for log2fc calculation
            cat > "${OUTPUT_DIR}/log2fc/calculate_log2fc_${treat_dox}_vs_${dam_dox}.R" << EOF
# Log2FC calculation for ${treat_dox} vs ${dam_dox}
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)

# Load binned data
treat_bins <- read.table("${OUTPUT_DIR}/analysis/${treat_dox}_bin_counts.bed", 
                         col.names=c("chr", "start", "end", "count"))
control_bins <- read.table("${OUTPUT_DIR}/analysis/${dam_dox}_bin_counts.bed", 
                           col.names=c("chr", "start", "end", "count"))

# Add pseudocount and calculate log2FC
treat_bins\$count_norm <- treat_bins\$count + 1
control_bins\$count_norm <- control_bins\$count + 1

# Merge data
merged_data <- data.frame(
  chr = treat_bins\$chr,
  start = treat_bins\$start,
  end = treat_bins\$end,
  treat_count = treat_bins\$count,
  control_count = control_bins\$count,
  log2FC = log2(treat_bins\$count_norm / control_bins\$count_norm)
)

# Write output
write.table(merged_data, file="${OUTPUT_DIR}/log2fc/${treat_dox}_vs_${dam_dox}_log2fc.txt",
            quote=FALSE, sep="\t", row.names=FALSE)

# Create bedGraph for visualization
bedgraph_data <- data.frame(
  chr = merged_data\$chr,
  start = merged_data\$start,
  end = merged_data\$end,
  score = merged_data\$log2FC
)

# Filter out extreme values
bedgraph_data\$score[bedgraph_data\$score > 5] <- 5
bedgraph_data\$score[bedgraph_data\$score < -5] <- -5

# Write bedGraph
write.table(bedgraph_data, file="${OUTPUT_DIR}/log2fc/${treat_dox}_vs_${dam_dox}_log2fc.bedgraph",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Create BigWig if possible
if (require(rtracklayer)) {
  try({
    gr <- GRanges(
      seqnames = bedgraph_data\$chr,
      ranges = IRanges(start = bedgraph_data\$start + 1, end = bedgraph_data\$end),
      score = bedgraph_data\$score
    )
    export(gr, "${OUTPUT_DIR}/log2fc/${treat_dox}_vs_${dam_dox}_log2fc.bw")
  })
}

cat("Log2FC calculation complete for ${treat_dox} vs ${dam_dox}\\n")
EOF
            
            # Run R script
            Rscript "${OUTPUT_DIR}/log2fc/calculate_log2fc_${treat_dox}_vs_${dam_dox}.R"
        done
    done
    
    log "Log2 fold change calculations complete."
}

create_sample_file() {
    log "No sample file provided. Creating template at ${SAMPLE_FILE}..."
    
    cat > "${SAMPLE_FILE}" << EOF
DAM-1
DAM-2
DAM-3
DAM-Dox-1
DAM-Dox-2
DAM-Dox-3
Gli3-Hand2-1
Gli3-Hand2-2
Gli3-Hand2-3
Gli3-Hand2-Dox-1
Gli3-Hand2-Dox-2
Gli3-Hand2-Dox-3
EOF
    
    log "Sample file template created. Edit this file with your sample names."
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
            shift 2
            ;;
        --fasta)
            GENOME_FASTA="$2"
            shift 2
            ;;
        --bowtie2)
            BOWTIE2_INDEX="$2"
            shift 2
            ;;
        --gatc)
            GATC_BED="$2"
            shift 2
            ;;
        --blacklist)
            BLACKLIST="$2"
            shift 2
            ;;
        --extension-script)
            EXTENSION_SCRIPT="$2"
            shift 2
            ;;
        --bin-size)
            BIN_SIZE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --genome-size)
            GENOME_SIZE="$2"
            shift 2
            ;;
        --control-pattern)
            DAM_ONLY_PATTERN="$2"
            shift 2
            ;;
        --treat-pattern)
            TREAT_PATTERN="$2"
            shift 2
            ;;
        --dox-pattern)
            DOX_PATTERN="$2"
            shift 2
            ;;
        --skip-log2fc)
            CALCULATE_LOG2FC=false
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

# Check if sample file exists, create template if not
if [ ! -f "${SAMPLE_FILE}" ]; then
    create_sample_file
    log "Please edit the sample file and run the script again."
    exit 1
fi

# Check dependencies
check_dependencies

# Create output directories
create_directories

# Create genome bins
create_genome_bins

# Process each sample
while read sample; do
    process_sample "${sample}"
done < "${SAMPLE_FILE}"

# Calculate log2 fold changes if requested
if [ "${CALCULATE_LOG2FC}" = true ]; then
    calculate_log2fc
fi

log "Split-DamID analysis workflow completed successfully!"
log "Results can be found in: ${OUTPUT_DIR}"

# Print summary
echo ""
echo "Summary of results:"
echo "==================="
echo "Processed $(wc -l < ${SAMPLE_FILE}) samples"
echo "BAM files: ${OUTPUT_DIR}/bam/"
echo "BigWig files: ${OUTPUT_DIR}/bigwig/"
echo "Binned count data: ${OUTPUT_DIR}/analysis/"
echo "Peak calls: ${OUTPUT_DIR}/peaks/"

if [ "${CALCULATE_LOG2FC}" = true ]; then
    echo "Log2 fold change data: ${OUTPUT_DIR}/log2fc/"
fi

echo ""
log "Done."
exit 0
