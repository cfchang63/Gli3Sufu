#!/usr/bin/env bash
#
# Split-DamID Complete Analysis Workflow
# 
# This script automates the entire Split-DamID analysis pipeline from
# raw data processing through binning, differential analysis, and peak annotation.
#
# Author: Ching-Fang Chang
# Date: March 23, 2025
#

set -e  # Exit on error
set -o pipefail  # Exit if any command in pipe fails

# CONFIGURATION
##############################################################################

# Default parameters (can be overridden via command line arguments)
RAW_DATA_DIR="./raw_data"
OUTPUT_DIR="./splitdamid_results"
FINAL_OUTPUT_DIR="./splitdamid_analysis"
GENOME="mm10"
GENOME_FASTA="${HOME}/Databases/${GENOME}/${GENOME}.fa"
BOWTIE2_INDEX="${HOME}/Databases/${GENOME}/${GENOME}_bt2index"
GATC_BED="${HOME}/Databases/${GENOME}/${GENOME}_GATC_sites.bed"
BLACKLIST="${HOME}/Databases/${GENOME}/blacklist_mtDNA.bed"
SAMPLE_FILE="samples.txt"
BIN_SIZE=75
THREADS=12
CONTRAST=""

# Paths to scripts
SPLITDAMID_SCRIPT="./splitdamid-workflow.sh"
BINNING_SCRIPT="./prepare_splitdamid_bins.py"
DIFFBIND_SCRIPT="./splitdamid_diff_analysis.R"

# FUNCTIONS
##############################################################################

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  --raw-data DIR            Directory with raw FASTQ files (default: ${RAW_DATA_DIR})"
    echo "  -o, --output DIR          Main output directory (default: ${OUTPUT_DIR})"
    echo "  --final-output DIR        Directory for final analysis results (default: ${FINAL_OUTPUT_DIR})"
    echo "  -g, --genome NAME         Genome name (default: ${GENOME})"
    echo "  --fasta FILE              Genome FASTA file (default: ${GENOME_FASTA})"
    echo "  --bowtie2 DIR             Bowtie2 index directory (default: ${BOWTIE2_INDEX})"
    echo "  --gatc FILE               GATC sites BED file (default: ${GATC_BED})"
    echo "  --blacklist FILE          Blacklist+mtDNA regions BED file (default: ${BLACKLIST})"
    echo "  -s, --samples FILE        Text file with sample names (default: ${SAMPLE_FILE})"
    echo "  --bin-size NUM            Size of genome bins in bp (default: ${BIN_SIZE})"
    echo "  -t, --threads NUM         Number of threads (default: ${THREADS})"
    echo "  --splitdamid-script FILE  Path to Split-DamID workflow script (default: ${SPLITDAMID_SCRIPT})"
    echo "  --binning-script FILE     Path to binning script (default: ${BINNING_SCRIPT})"
    echo "  --diffbind-script FILE    Path to differential binding analysis script (default: ${DIFFBIND_SCRIPT})"
    echo "  --contrast NUM            Specific contrast to analyze (optional)"
    echo "  --only-process            Only process raw data (skip differential analysis)"
    echo "  --only-analyze            Skip raw data processing, only do differential analysis"
    echo "  -h, --help                Display this help message"
    exit 1
}

log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

check_dependencies() {
    local dependencies=(python3 Rscript bowtie2 samtools bedtools)
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
    
    # Check script files
    if [ ! -f "${SPLITDAMID_SCRIPT}" ]; then
        log "ERROR: Split-DamID workflow script not found at ${SPLITDAMID_SCRIPT}"
        log "Please specify the correct path with --splitdamid-script"
        exit 1
    fi
    
    if [ ! -f "${BINNING_SCRIPT}" ]; then
        log "ERROR: Binning script not found at ${BINNING_SCRIPT}"
        log "Please specify the correct path with --binning-script"
        exit 1
    fi
    
    if [ ! -f "${DIFFBIND_SCRIPT}" ]; then
        log "ERROR: Differential binding analysis script not found at ${DIFFBIND_SCRIPT}"
        log "Please specify the correct path with --diffbind-script"
        exit 1
    fi
    
    # Set script permissions
    chmod +x "${SPLITDAMID_SCRIPT}"
    chmod +x "${DIFFBIND_SCRIPT}"
    
    log "All dependencies found."
}

create_directories() {
    mkdir -p "${OUTPUT_DIR}"
    mkdir -p "${FINAL_OUTPUT_DIR}"
    log "Created output directories: ${OUTPUT_DIR}, ${FINAL_OUTPUT_DIR}"
}

process_raw_data() {
    log "Processing raw Split-DamID data..."
    
    # Check if sample file exists
    if [ ! -f "${SAMPLE_FILE}" ]; then
        log "ERROR: Sample file not found: ${SAMPLE_FILE}"
        exit 1
    fi
    
    # Run Split-DamID workflow
    "${SPLITDAMID_SCRIPT}" \
        --samples "${SAMPLE_FILE}" \
        --output "${OUTPUT_DIR}" \
        --genome "${GENOME}" \
        --fasta "${GENOME_FASTA}" \
        --bowtie2 "${BOWTIE2_INDEX}" \
        --gatc "${GATC_BED}" \
        --blacklist "${BLACKLIST}" \
        --bin-size "${BIN_SIZE}" \
        --threads "${THREADS}"
    
    log "Raw data processing complete."
}

prepare_binned_counts() {
    log "Preparing binned count matrix..."
    
    # Check if bin counts exist
    bin_files=$(find "${OUTPUT_DIR}/analysis" -name "*_bin_counts.bed" | wc -l)
    
    if [ "${bin_files}" -eq 0 ]; then
        log "ERROR: No bin count files found in ${OUTPUT_DIR}/analysis"
        exit 1
    fi
    
    # Create sample info file
    if [ ! -d "${FINAL_OUTPUT_DIR}" ]; then
        mkdir -p "${FINAL_OUTPUT_DIR}"
    fi
    
    # Run binning script
    python3 "${BINNING_SCRIPT}" \
        --input "${OUTPUT_DIR}" \
        --output "${FINAL_OUTPUT_DIR}/binned_count_matrix.csv" \
        --sample-info "${FINAL_OUTPUT_DIR}/sample_info.csv" \
        --bin-dir analysis \
        --bin-pattern "*_bin_counts.bed"
    
    log "Binned count matrix created: ${FINAL_OUTPUT_DIR}/binned_count_matrix.csv"
    log "Sample info file: ${FINAL_OUTPUT_DIR}/sample_info.csv"
}

run_differential_analysis() {
    log "Running differential binding analysis..."
    
    # Check if count matrix and sample info exist
    if [ ! -f "${FINAL_OUTPUT_DIR}/binned_count_matrix.csv" ]; then
        log "ERROR: Count matrix not found: ${FINAL_OUTPUT_DIR}/binned_count_matrix.csv"
        exit 1
    fi
    
    if [ ! -f "${FINAL_OUTPUT_DIR}/sample_info.csv" ]; then
        log "ERROR: Sample info file not found: ${FINAL_OUTPUT_DIR}/sample_info.csv"
        exit 1
    fi
    
    # Prepare contrast argument
    contrast_arg=""
    if [ -n "${CONTRAST}" ]; then
        contrast_arg="${CONTRAST}"
    fi
    
    # Run differential analysis script
    Rscript "${DIFFBIND_SCRIPT}" \
        "${FINAL_OUTPUT_DIR}/binned_count_matrix.csv" \
        "${FINAL_OUTPUT_DIR}/sample_info.csv" \
        "${FINAL_OUTPUT_DIR}" \
        ${contrast_arg}
    
    log "Differential binding analysis complete. Results in ${FINAL_OUTPUT_DIR}"
}

create_integrated_report() {
    log "Creating integrated analysis report..."
    
    # Create report directory
    mkdir -p "${FINAL_OUTPUT_DIR}/report"
    
    # Create HTML report
    cat > "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html" << EOF
<!DOCTYPE html>
<html>
<head>
    <title>Split-DamID Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; line-height: 1.6; margin: 40px; }
        h1, h2, h3 { color: #2c3e50; }
        .container { max-width: 1200px; margin: 0 auto; }
        .section { margin-bottom: 30px; border-bottom: 1px solid #eee; padding-bottom: 20px; }
        .figure { margin: 20px 0; text-align: center; }
        .figure img { max-width: 100%; border: 1px solid #ddd; }
        .caption { font-style: italic; color: #666; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #f2f2f2; }
        tr:hover { background-color: #f5f5f5; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Split-DamID Analysis Report</h1>
        <p>Generated on: $(date)</p>
        
        <div class="section">
            <h2>Sample Information</h2>
            <p>The following samples were analyzed in this Split-DamID experiment:</p>
            <table>
                <tr>
                    <th>Sample</th>
                    <th>Group</th>
                    <th>Replicate</th>
                </tr>
EOF

    # Add sample information to the report
    if [ -f "${FINAL_OUTPUT_DIR}/sample_info.csv" ]; then
        # Skip header line
        tail -n +2 "${FINAL_OUTPUT_DIR}/sample_info.csv" | while IFS=, read -r sample group replicate; do
            echo "<tr><td>${sample}</td><td>${group}</td><td>${replicate}</td></tr>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
        done
    fi

    # Continue HTML report
    cat >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html" << EOF
            </table>
        </div>
        
        <div class="section">
            <h2>Quality Control</h2>
            
            <div class="figure">
                <h3>Principal Component Analysis</h3>
                <img src="../PCA_plot.pdf" alt="PCA Plot">
                <p class="caption">Principal Component Analysis showing sample clustering.</p>
            </div>
            
            <div class="figure">
                <h3>Sample Correlation</h3>
                <img src="../sample_distance_heatmap.pdf" alt="Sample Distance Heatmap">
                <p class="caption">Heatmap showing sample-to-sample distances based on binding profiles.</p>
            </div>
        </div>
        
        <div class="section">
            <h2>Differential Binding Analysis</h2>
            <p>Summary of differential binding results:</p>
            
            <table>
                <tr>
                    <th>Comparison</th>
                    <th>Total Bins</th>
                    <th>Significant Bins (FDR < 0.05)</th>
                    <th>Upregulated</th>
                    <th>Downregulated</th>
                </tr>
EOF

    # Add differential binding results to the report
    for results_file in "${FINAL_OUTPUT_DIR}"/significant_bins_*.csv; do
        if [ -f "${results_file}" ]; then
            # Extract comparison name from filename
            comparison=$(basename "${results_file}" | sed 's/significant_bins_\(.*\)\.csv/\1/')
            
            # Count total, up, and down regulated bins
            total_bins=$(wc -l < "${results_file}")
            total_bins=$((total_bins - 1))  # Subtract header line
            
            # Count upregulated bins (log2FoldChange > 0)
            up_bins=$(awk -F, 'NR>1 && $3 > 0 {count++} END {print count}' "${results_file}")
            if [ -z "${up_bins}" ]; then up_bins=0; fi
            
            # Count downregulated bins (log2FoldChange < 0)
            down_bins=$(awk -F, 'NR>1 && $3 < 0 {count++} END {print count}' "${results_file}")
            if [ -z "${down_bins}" ]; then down_bins=0; fi
            
            echo "<tr><td>${comparison}</td><td>${total_bins}</td><td>${total_bins}</td><td>${up_bins}</td><td>${down_bins}</td></tr>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
        fi
    done

    # Continue HTML report
    cat >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html" << EOF
            </table>
            
            <h3>Volcano Plots</h3>
            <p>Volcano plots showing differential binding for each comparison:</p>
EOF

    # Add volcano plots to the report
    for volcano_file in "${FINAL_OUTPUT_DIR}"/volcano_plot_*.pdf; do
        if [ -f "${volcano_file}" ]; then
            # Extract comparison name from filename
            comparison=$(basename "${volcano_file}" | sed 's/volcano_plot_\(.*\)\.pdf/\1/')
            
            echo "<div class=\"figure\">" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            echo "<h4>${comparison}</h4>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            echo "<img src=\"../${comparison}/volcano_plot_${comparison}.pdf\" alt=\"Volcano Plot: ${comparison}\">" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            echo "<p class=\"caption\">Volcano plot showing differentially bound regions for ${comparison}.</p>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            echo "</div>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
        fi
    done

    # Continue HTML report
    cat >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html" << EOF
        </div>
        
        <div class="section">
            <h2>Functional Enrichment Analysis</h2>
            <p>Gene Ontology and KEGG pathway enrichment results:</p>
EOF

    # Add GO and KEGG plots to the report
    for go_file in "${FINAL_OUTPUT_DIR}"/GO_dotplot_*.pdf; do
        if [ -f "${go_file}" ]; then
            # Extract comparison name from filename
            comparison=$(basename "${go_file}" | sed 's/GO_dotplot_\(.*\)\.pdf/\1/')
            
            echo "<h3>Enrichment for ${comparison}</h3>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            
            echo "<div class=\"figure\">" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            echo "<h4>Gene Ontology Enrichment</h4>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            echo "<img src=\"../GO_dotplot_${comparison}.pdf\" alt=\"GO Enrichment: ${comparison}\">" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            echo "<p class=\"caption\">Gene Ontology terms enriched in differential binding regions for ${comparison}.</p>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            echo "</div>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            
            if [ -f "${FINAL_OUTPUT_DIR}/KEGG_dotplot_${comparison}.pdf" ]; then
                echo "<div class=\"figure\">" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
                echo "<h4>KEGG Pathway Enrichment</h4>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
                echo "<img src=\"../KEGG_dotplot_${comparison}.pdf\" alt=\"KEGG Enrichment: ${comparison}\">" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
                echo "<p class=\"caption\">KEGG pathways enriched in differential binding regions for ${comparison}.</p>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
                echo "</div>" >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
            fi
        fi
    done

    # Finish HTML report
    cat >> "${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html" << EOF
        </div>
        
        <div class="section">
            <h2>Conclusion</h2>
            <p>This report summarizes the results of a Split-DamID analysis workflow, which identified differentially bound regions between the compared conditions.</p>
            <p>The full results, including detailed tables of differential binding and functional enrichment, can be found in the analysis directory.</p>
        </div>
        
        <div class="section">
            <h2>Methods</h2>
            <p>The analysis was performed using the following workflow:</p>
            <ol>
                <li>Raw Split-DamID data was processed with the Split-DamID workflow script</li>
                <li>Reads were aligned to the ${GENOME} genome</li>
                <li>Aligned reads were extended to the nearest GATC sites</li>
                <li>The genome was binned into ${BIN_SIZE}bp bins</li>
                <li>Reads were counted in each bin for each sample</li>
                <li>Differential binding analysis was performed using DESeq2</li>
                <li>Differential bins were annotated to nearby genes</li>
                <li>Functional enrichment analysis was performed for annotated genes</li>
            </ol>
        </div>
        
        <div>
            <p>Report generated by Split-DamID Complete Analysis Workflow</p>
            <p>Date: $(date)</p>
        </div>
    </div>
</body>
</html>
EOF

    log "Integrated analysis report created: ${FINAL_OUTPUT_DIR}/report/splitdamid_analysis_report.html"
}

# MAIN SCRIPT
##############################################################################

# Variables to control workflow steps
ONLY_PROCESS=false
ONLY_ANALYZE=false

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --raw-data)
            RAW_DATA_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --final-output)
            FINAL_OUTPUT_DIR="$2"
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
        -s|--samples)
            SAMPLE_FILE="$2"
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
        --splitdamid-script)
            SPLITDAMID_SCRIPT="$2"
            shift 2
            ;;
        --binning-script)
            BINNING_SCRIPT="$2"
            shift 2
            ;;
        --diffbind-script)
            DIFFBIND_SCRIPT="$2"
            shift 2
            ;;
        --contrast)
            CONTRAST="$2"
            shift 2
            ;;
        --only-process)
            ONLY_PROCESS=true
            shift
            ;;
        --only-analyze)
            ONLY_ANALYZE=true
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

# Check for conflicting options
if [ "${ONLY_PROCESS}" = true ] && [ "${ONLY_ANALYZE}" = true ]; then
    log "ERROR: Cannot specify both --only-process and --only-analyze"
    exit 1
fi

# Check dependencies
check_dependencies

# Create output directories
create_directories

# Main workflow steps
if [ "${ONLY_ANALYZE}" = false ]; then
    # Process raw data
    process_raw_data
fi

if [ "${ONLY_PROCESS}" = false ]; then
    # Prepare binned counts
    prepare_binned_counts
    
    # Run differential analysis
    run_differential_analysis
    
    # Create integrated report
    create_integrated_report
fi

log "Split-DamID complete analysis workflow completed successfully!"
log "Results can be found in: ${FINAL_OUTPUT_DIR}"

exit 0