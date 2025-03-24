# Split-DamID Differential Binding Analysis

This document explains how to use the `splitdamid_diffbind_analysis.R` script to perform differential binding analysis of Split-DamID data.

## Overview

After processing raw Split-DamID data through alignment, GATC extension, and peak calling, this script performs:

1. **Differential Binding Analysis**: Identifies differentially bound regions between conditions
2. **Peak Annotation**: Annotates peaks with genomic features (promoters, introns, etc.)
3. **Functional Enrichment**: Performs GO and KEGG pathway enrichment analysis

## Prerequisites

The following R packages are required:
- DiffBind
- ChIPseeker
- clusterProfiler
- TxDb.Mmusculus.UCSC.mm10.knownGene
- org.Mm.eg.db
- ggupset
- EnhancedVolcano

The script will check for these packages and install them if they're missing.

## Sample Sheet Format

The script requires a DiffBind sample sheet in CSV format with the following columns:

```
SampleID,Tissue,Factor,Condition,Replicate,bamReads,Peaks,PeakCaller
DAM-1,NIH3T3,DAM,DAM,1,/path/to/DAM-1_filtered_sorted.bam,/path/to/DAM-1_broad_peaks.broadPeak,bed
DAM-2,NIH3T3,DAM,DAM,2,/path/to/DAM-2_filtered_sorted.bam,/path/to/DAM-2_broad_peaks.broadPeak,bed
...
```

Key columns:
- **SampleID**: Unique identifier for each sample
- **Condition**: Group identifier (e.g., DAM, Gli3-Hand2, DAM-Dox, Gli3-Hand2-Dox)
- **Replicate**: Replicate number within the condition
- **bamReads**: Path to filtered BAM file
- **Peaks**: Path to peak file from MACS3
- **PeakCaller**: Set to "bed" for MACS3 broad peaks

A template is provided in `DiffBindSampleSheet.csv`.

## Usage

```bash
Rscript splitdamid_diffbind_analysis.R <samplesheet.csv> [output_dir] [contrast_number]
```

### Arguments:

1. **samplesheet.csv** (required): Path to the DiffBind sample sheet
2. **output_dir** (optional): Directory for output files (default: `./diffbind_results`)
3. **contrast_number** (optional): Number of the contrast to analyze (default: 2)

### Example:

```bash
Rscript splitdamid_diffbind_analysis.R DiffBindSampleSheet.csv ./results 2
```

This will analyze the contrast number 2 (usually DAM vs Gli3-Hand2) and output results to the `./results` directory.

## Available Contrasts

The script shows all available contrasts at runtime. The default contrasts are:

1. DAM vs DAM-Dox
2. DAM vs Gli3-Hand2
3. DAM vs Gli3-Hand2-Dox
4. DAM-Dox vs Gli3-Hand2
5. DAM-Dox vs Gli3-Hand2-Dox
6. Gli3-Hand2-Dox vs Gli3-Hand2

You can specify which contrast to analyze using the third command line argument.

## Output Files

The script generates the following output:

### Main Results
- `<contrast>_diffbind_results.csv`: Table of differential peaks with statistics
- `<contrast>_annotated_peaks.csv`: Peaks annotated with genomic features and genes
- `<contrast>_diffbind_analysis.RData`: R workspace with full analysis results

### Plots Directory
- `correlation_heatmap.pdf`: Correlation between samples
- `<contrast>_MA_plot.pdf`: MA plot of differential binding
- `<contrast>_volcano_plot.pdf`: Volcano plot
- `<contrast>_differential_binding_heatmap.pdf`: Heatmap of differentially bound regions
- `<contrast>_PCA_plot.pdf`: Principal Component Analysis
- `<contrast>_peak_annotation_pie.pdf`: Pie chart of peak annotations
- `<contrast>_peak_annotation_bar.pdf`: Bar chart of peak annotations
- `<contrast>_peak_distance_to_TSS.pdf`: Distance to nearest TSS
- `<contrast>_enhanced_volcano_plot.pdf`: Enhanced volcano plot with gene labels
- `<contrast>_upset_plot.pdf`: UpSet plot showing overlapping features

### Enrichment Directory
- `<contrast>_GO_enrichment.csv`: Gene Ontology enrichment results
- `<contrast>_GO_enrichment_barplot.pdf`: Bar plot of GO terms
- `<contrast>_GO_enrichment_dotplot.pdf`: Dot plot of GO terms
- `<contrast>_GO_enrichment_map.pdf`: Enrichment map for GO terms
- `<contrast>_KEGG_enrichment.csv`: KEGG pathway enrichment results
- `<contrast>_KEGG_enrichment_barplot.pdf`: Bar plot of KEGG pathways
- `<contrast>_KEGG_enrichment_dotplot.pdf`: Dot plot of KEGG pathways
- `<contrast>_KEGG_enrichment_map.pdf`: Enrichment map for KEGG pathways

## Creating the Sample Sheet

After running the Split-DamID workflow, you can create the sample sheet with:

```bash
# Base path to your Split-DamID results
BASE_DIR="splitdamid_results"

# Create sample sheet
echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,Peaks,PeakCaller" > DiffBindSampleSheet.csv

# Add each sample
for SAMPLE in DAM-1 DAM-2 DAM-3 DAM-Dox-1 DAM-Dox-2 DAM-Dox-3 Gli3-Hand2-1 Gli3-Hand2-2 Gli3-Hand2-3 Gli3-Hand2-Dox-1 Gli3-Hand2-Dox-2 Gli3-Hand2-Dox-3; do
  # Extract condition and replicate
  if [[ $SAMPLE == *"Gli3-Hand2-Dox"* ]]; then
    CONDITION="Gli3-Hand2-Dox"
    REPLICATE=${SAMPLE##*-}
    FACTOR="Gli3-Hand2"
  elif [[ $SAMPLE == *"Gli3-Hand2"* ]]; then
    CONDITION="Gli3-Hand2"
    REPLICATE=${SAMPLE##*-}
    FACTOR="Gli3-Hand2"
  elif [[ $SAMPLE == *"DAM-Dox"* ]]; then
    CONDITION="DAM-Dox"
    REPLICATE=${SAMPLE##*-}
    FACTOR="DAM"
  else
    CONDITION="DAM"
    REPLICATE=${SAMPLE##*-}
    FACTOR="DAM"
  fi
  
  echo "$SAMPLE,NIH3T3,$FACTOR,$CONDITION,$REPLICATE,$BASE_DIR/bam/${SAMPLE}_filtered_sorted.bam,$BASE_DIR/peaks/${SAMPLE}/${SAMPLE}_broad_peaks.broadPeak,bed" >> DiffBindSampleSheet.csv
done
```

## Interpreting Results

### Differential Binding Results
The `<contrast>_diffbind_results.csv` file contains:
- **Chromosome, Start, End**: Genomic coordinates of the peak
- **Fold**: Log2 fold change (positive = enriched in condition 1, negative = enriched in condition 2)
- **p-value**: Statistical significance of the difference
- **FDR**: False Discovery Rate (adjusted p-value)

### Peak Annotation
The `<contrast>_annotated_peaks.csv` file contains:
- Genomic annotation (promoter, intron, exon, etc.)
- Nearest gene information
- Distance to nearest TSS

### Functional Enrichment
The GO and KEGG enrichment files contain:
- Enriched terms/pathways
- Gene counts
- P-values and q-values
- Gene ratios

## Troubleshooting

1. **Missing peaks or BAM files**:
   - The script will warn if any peak or BAM files are missing
   - Verify file paths in the sample sheet

2. **No significant peaks**:
   - Try adjusting the FDR threshold in the script (default: 0.05)
   - Check if the correct contrast is being analyzed

3. **R package installation issues**:
   - Make sure Bioconductor is properly installed
   - Run `BiocManager::install()` to update Bioconductor

4. **Memory issues**:
   - For large datasets, increase R's memory limit: `R --max-mem-size=16G`

## References

1. Ross-Innes, C.S., et al. (2012). Differential oestrogen receptor binding is associated with clinical outcome in breast cancer. Nature, 481, 389-393.
2. Yu, G., et al. (2015). ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics, 31, 2382-2383.
3. Yu, G., et al. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16, 284-287.
