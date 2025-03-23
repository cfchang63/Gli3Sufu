#!/usr/bin/env Rscript

#' ATAC-seq Differential Analysis and Visualization
#' 
#' This script performs differential accessibility analysis on ATAC-seq data
#' and creates various visualizations to help interpret the results.
#' 
#' Usage:
#'   Rscript analyze_atac.R --counts counts_matrix.txt --metadata metadata.txt --output results_dir
#'
#' Author: [Your Name]
#' Date: March 23, 2025

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(DiffBind)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(clusterProfiler)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--counts"), type="character", default=NULL,
              help="Path to counts matrix file"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Path to metadata file"),
  make_option(c("--output"), type="character", default="atac_results",
              help="Output directory [default=%default]"),
  make_option(c("--fdr"), type="numeric", default=0.05,
              help="FDR threshold for differential peaks [default=%default]"),
  make_option(c("--log2fc"), type="numeric", default=1.0,
              help="Log2 fold change threshold [default=%default]"),
  make_option(c("--genome"), type="character", default="mm10",
              help="Genome assembly [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$counts) || is.null(opt$metadata)) {
  stop("Both --counts and --metadata arguments are required")
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Set up logging
log_file <- file.path(opt$output, "analysis_log.txt")
log_conn <- file(log_file, open = "w")
log <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_message <- paste0("[", timestamp, "] ", message)
  cat(full_message, "\n", file = log_conn, append = TRUE)
  cat(full_message, "\n")
}

# Log script start
log("Starting ATAC-seq differential analysis")
log(paste("Counts file:", opt$counts))
log(paste("Metadata file:", opt$metadata))
log(paste("Output directory:", opt$output))

# Load and prepare data
#--------------------------------------------------

# Load metadata
log("Loading metadata")
metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE)
log(paste("Found", nrow(metadata), "samples in metadata"))

# Prepare DiffBind sample sheet
create_sample_sheet <- function(metadata) {
  # This assumes standard file naming from the ATAC-seq pipeline
  sample_sheet <- data.frame(
    SampleID = metadata$sample_id,
    Condition = metadata$condition,
    Replicate = metadata$replicate,
    bamReads = file.path("filtered", paste0(metadata$sample_id, ".shifted.sorted.bam")),
    Peaks = file.path("peaks", paste0(metadata$sample_id), paste0(metadata$sample_id, "_peaks.narrowPeak")),
    PeakCaller = "narrowPeak",
    stringsAsFactors = FALSE
  )
  
  # Add control samples if available
  if ("control" %in% colnames(metadata) && any(!is.na(metadata$control))) {
    sample_sheet$ControlID <- metadata$control
    sample_sheet$bamControl <- file.path("filtered", paste0(metadata$control, ".shifted.sorted.bam"))
  }
  
  return(sample_sheet)
}

# Create sample sheet
sample_sheet <- create_sample_sheet(metadata)
sample_sheet_file <- file.path(opt$output, "sample_sheet.csv")
write.csv(sample_sheet, sample_sheet_file, row.names = FALSE)
log(paste("Created DiffBind sample sheet:", sample_sheet_file))

# Load count matrix (if provided) or create from sample sheet
if (!is.null(opt$counts) && file.exists(opt$counts)) {
  log("Loading pre-computed count matrix")
  counts_matrix <- read.delim(opt$counts, row.names = 1)
  
  # Create DiffBind object from counts
  dba <- dba.peakset(NULL, sampleSheet = sample_sheet)
  dba_counts <- dba.count(dba, score = DBA_SCORE_READS, 
                          peaks = NULL, # Use all peaks in sample sheet
                          readFormat = "bam", # Assuming BAM format
                          summits = FALSE)
} else {
  log("Creating DiffBind object and counting reads")
  dba <- dba.peakset(sampleSheet = sample_sheet, 
                    dir = getwd(), # Assumes paths in sample sheet are relative to working dir
                    minOverlap = 2) # Require peaks to be in at least 2 samples
  
  log(paste("Found", dba$config$peakCaller, "peaks"))
  log(paste("Total peak count:", length(dba$peaks[[1]])))
  
  # Count reads in peaks
  log("Counting reads in consensus peaks")
  dba_counts <- dba.count(dba, score = DBA_SCORE_READS, 
                          summits = FALSE, # Don't re-center on summits
                          filter = 1) # Filter out peaks with no reads
  
  # Export count matrix
  counts_matrix <- dba.peakset(dba_counts, bRetrieve = TRUE)
  counts_file <- file.path(opt$output, "peak_counts.txt")
  write.table(counts_matrix, counts_file, sep = "\t", quote = FALSE)
  log(paste("Saved count matrix to:", counts_file))
}

# Normalize and explore data
#--------------------------------------------------

# Normalize data
log("Performing normalization")
dba_norm <- dba.normalize(dba_counts, normalize = DBA_NORM_LIB)

# Generate PCA plot
log("Generating PCA plot")
pdf(file.path(opt$output, "PCA_plot.pdf"), width = 8, height = 6)
dba.plotPCA(dba_norm, attributes = DBA_CONDITION, label = DBA_ID)
dev.off()

# Generate correlation heatmap
log("Generating correlation heatmap")
pdf(file.path(opt$output, "correlation_heatmap.pdf"), width = 8, height = 7)
dba.plotHeatmap(dba_norm, correlations = TRUE, 
                colSideCols = c("blue", "red")[as.numeric(factor(dba_norm$samples$Condition))],
                cexCol = 0.8, cexRow = 0.8)
dev.off()

# Differential analysis
#--------------------------------------------------

# Get all condition pairs for comparison
conditions <- unique(metadata$condition)
log(paste("Found conditions:", paste(conditions, collapse = ", ")))

comparisons <- combn(conditions, 2, simplify = FALSE)
log(paste("Will perform", length(comparisons), "pairwise comparisons"))

# Perform differential analysis for each comparison
for (comp in comparisons) {
  cond1 <- comp[1]
  cond2 <- comp[2]
  comparison_name <- paste(cond1, "vs", cond2, sep = "_")
  
  log(paste("Performing differential analysis:", comparison_name))
  
  # Create contrast
  dba_contrast <- dba.contrast(dba_norm, 
                              group1 = dba_norm$samples$Condition %in% cond1,
                              group2 = dba_norm$samples$Condition %in% cond2,
                              name1 = cond1, 
                              name2 = cond2)
  
  # Analyze contrast
  dba_analyze <- dba.analyze(dba_contrast, method = DBA_DESEQ2)
  
  # Retrieve results
  dba_results <- dba.report(dba_analyze, 
                           contrast = 1, 
                           th = 1) # Return all results, filter later
  
  # Create results directory for this comparison
  comp_dir <- file.path(opt$output, comparison_name)
  if (!dir.exists(comp_dir)) dir.create(comp_dir)
  
  # Export results table
  results_df <- as.data.frame(dba_results)
  results_file <- file.path(comp_dir, paste0(comparison_name, "_diffpeaks.tsv"))
  write.table(results_df, results_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Filter significant peaks
  sig_peaks <- results_df %>%
    filter(FDR < opt$fdr, abs(Fold) > opt$log2fc)
  
  log(paste("Found", nrow(sig_peaks), "significant differential peaks (FDR <", 
            opt$fdr, ", |log2FC| >", opt$log2fc, ")"))
  
  # Create MA plot
  pdf(file.path(comp_dir, paste0(comparison_name, "_MA_plot.pdf")), width = 7, height = 6)
  dba.plotMA(dba_analyze, contrast = 1, th = opt$fdr, bXY = TRUE)
  dev.off()
  
  # Create volcano plot
  pdf(file.path(comp_dir, paste0(comparison_name, "_volcano_plot.pdf")), width = 7, height = 6)
  
  # Custom volcano plot with ggplot2
  ggplot(results_df, aes(x = Fold, y = -log10(FDR))) +
    geom_point(aes(color = ifelse(FDR < opt$fdr & abs(Fold) > opt$log2fc, 
                                 ifelse(Fold > 0, "Up", "Down"), "NS")), 
               alpha = 0.7, size = 1) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey50")) +
    geom_hline(yintercept = -log10(opt$fdr), linetype = "dashed") +
    geom_vline(xintercept = c(-opt$log2fc, opt$log2fc), linetype = "dashed") +
    labs(x = "Log2 Fold Change", y = "-Log10 FDR", 
         title = paste("Differential Accessibility:", comparison_name),
         color = "Regulation") +
    theme_bw() +
    theme(legend.position = "top")
  
  dev.off()
  
  # Create BED files for each direction
  if (nrow(sig_peaks) > 0) {
    # Up-regulated peaks
    up_peaks <- sig_peaks %>% filter(Fold > 0)
    if (nrow(up_peaks) > 0) {
      up_bed <- up_peaks %>% 
        select(seqnames, start, end, width, Fold, FDR) %>%
        mutate(name = paste0("peak_", row_number()),
               score = -log10(FDR) * 10,
               strand = ".") %>%
        select(seqnames, start, end, name, score, strand, Fold, FDR)
      
      write.table(up_bed, file.path(comp_dir, paste0(comparison_name, "_up.bed")),
                  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    
    # Down-regulated peaks
    down_peaks <- sig_peaks %>% filter(Fold < 0)
    if (nrow(down_peaks) > 0) {
      down_bed <- down_peaks %>% 
        select(seqnames, start, end, width, Fold, FDR) %>%
        mutate(name = paste0("peak_", row_number()),
               score = -log10(FDR) * 10,
               strand = ".") %>%
        select(seqnames, start, end, name, score, strand, Fold, FDR)
      
      write.table(down_bed, file.path(comp_dir, paste0(comparison_name, "_down.bed")),
                  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  
  # Annotate peaks with nearest genes
  if (nrow(sig_peaks) > 0 && opt$genome == "mm10") {
    log("Annotating peaks with genomic features")
    
    # Convert to GRanges
    sig_gr <- makeGRangesFromDataFrame(sig_peaks, keep.extra.columns = TRUE)
    
    # Load TxDb object based on genome
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    
    # Annotate peaks
    peakAnno <- annotatePeak(sig_gr, TxDb = txdb,
                            annoDb = "org.Mm.eg.db",
                            tssRegion = c(-2000, 500))
    
    # Save annotation results
    anno_df <- as.data.frame(peakAnno)
    anno_file <- file.path(comp_dir, paste0(comparison_name, "_annotated.tsv"))
    write.table(anno_df, anno_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Plot annotation statistics
    pdf(file.path(comp_dir, paste0(comparison_name, "_annotation_plots.pdf")), width = 10, height = 8)
    plotAnnoPie(peakAnno)
    plotDistToTSS(peakAnno)
    dev.off()
    
    # Separate up and down regulated peaks for GO analysis
    if (sum(sig_peaks$Fold > 0) > 5) {
      up_genes <- anno_df %>% 
        filter(Fold > 0) %>% 
        pull(SYMBOL) %>% 
        unique()
      
      # Perform GO enrichment analysis
      log(paste("Performing GO enrichment analysis for", length(up_genes), "up-regulated genes"))
      ego_up <- enrichGO(gene = up_genes,
                         OrgDb = org.Mm.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1)
      
      if (!is.null(ego_up) && nrow(ego_up) > 0) {
        # Save GO results
        go_up_file <- file.path(comp_dir, paste0(comparison_name, "_up_GO.tsv"))
        write.table(as.data.frame(ego_up), go_up_file, sep = "\t", row.names = FALSE, quote = FALSE)
        
        # Plot GO results
        pdf(file.path(comp_dir, paste0(comparison_name, "_up_GO_dotplot.pdf")), width = 9, height = 7)
        dotplot(ego_up, showCategory = 15, title = paste("GO enrichment for up-regulated peaks in", comparison_name))
        dev.off()
      }
    }
    
    if (sum(sig_peaks$Fold < 0) > 5) {
      down_genes <- anno_df %>% 
        filter(Fold < 0) %>% 
        pull(SYMBOL) %>% 
        unique()
      
      log(paste("Performing GO enrichment analysis for", length(down_genes), "down-regulated genes"))
      ego_down <- enrichGO(gene = down_genes,
                           OrgDb = org.Mm.eg.db,
                           keyType = "SYMBOL",
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.1)
      
      if (!is.null(ego_down) && nrow(ego_down) > 0) {
        # Save GO results
        go_down_file <- file.path(comp_dir, paste0(comparison_name, "_down_GO.tsv"))
        write.table(as.data.frame(ego_down), go_down_file, sep = "\t", row.names = FALSE, quote = FALSE)
        
        # Plot GO results
        pdf(file.path(comp_dir, paste0(comparison_name, "_down_GO_dotplot.pdf")), width = 9, height = 7)
        dotplot(ego_down, showCategory = 15, title = paste("GO enrichment for down-regulated peaks in", comparison_name))
        dev.off()
      }
    }
  }
}

# Generate consolidated report
#--------------------------------------------------

# Create summary table of all comparisons
log("Generating summary report")

summary_data <- data.frame()
for (comp in comparisons) {
  cond1 <- comp[1]
  cond2 <- comp[2]
  comparison_name <- paste(cond1, "vs", cond2, sep = "_")
  
  # Check if results file exists
  results_file <- file.path(opt$output, comparison_name, paste0(comparison_name, "_diffpeaks.tsv"))
  if (file.exists(results_file)) {
    results <- read.delim(results_file)
    
    # Count significant peaks
    sig_up <- sum(results$FDR < opt$fdr & results$Fold > opt$log2fc)
    sig_down <- sum(results$FDR < opt$fdr & results$Fold < -opt$log2fc)
    
    # Add to summary table
    comp_summary <- data.frame(
      Comparison = comparison_name,
      Condition1 = cond1,
      Condition2 = cond2,
      Total_Peaks = nrow(results),
      Up_Peaks = sig_up,
      Down_Peaks = sig_down,
      Significant_Peaks = sig_up + sig_down,
      Percent_Significant = round((sig_up + sig_down) / nrow(results) * 100, 2)
    )
    
    summary_data <- rbind(summary_data, comp_summary)
  }
}

# Save summary table
if (nrow(summary_data) > 0) {
  summary_file <- file.path(opt$output, "differential_analysis_summary.tsv")
  write.table(summary_data, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
  log(paste("Saved summary report to:", summary_file))
}

# Close log file
close(log_conn)
log("ATAC-seq differential analysis complete!")
