#!/usr/bin/env Rscript

#' Split-DamID Differential Analysis and Peak Annotation
#' 
#' This script performs differential binding analysis on Split-DamID binned data,
#' followed by annotation of the differential bins to nearby genes.
#' 
#' Usage:
#'   Rscript splitdamid_diff_analysis.R <count_matrix> <sample_info> <output_dir> [contrast_number]
#'
#' Arguments:
#'   count_matrix    - Path to binned count matrix CSV file
#'   sample_info     - Path to sample information CSV file
#'   output_dir      - Directory for output files
#'   contrast_number - (Optional) Specific contrast to analyze (default: analyze all)
#'
#' Author: Ching-Fang Chang
#' Date: March 23, 2025

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(BiocParallel)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(ggrepel)
})

#' Check package installation and install if needed
check_and_install <- function(pkg_name, is_bioc = FALSE) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    message(paste0("Installing package: ", pkg_name))
    if (is_bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg_name)
    } else {
      install.packages(pkg_name)
    }
  }
}

#' Read and preprocess count matrix
read_count_matrix <- function(count_file) {
  message("Reading count matrix:", count_file)
  
  # Read the count matrix
  counts <- data.table::fread(count_file, data.table = FALSE)
  
  # Check if Bin_ID column exists
  if (!"Bin_ID" %in% colnames(counts)) {
    stop("Count matrix must contain a 'Bin_ID' column")
  }
  
  # Set Bin_ID as rownames and remove the column
  rownames(counts) <- counts$Bin_ID
  counts <- counts[, !colnames(counts) %in% c("Bin_ID", "Chromosome", "Start", "End")]
  
  # Check that we have count data
  if (ncol(counts) == 0) {
    stop("No count data found in count matrix after removing bin coordinates")
  }
  
  message(paste("Loaded count matrix with", nrow(counts), "bins and", ncol(counts), "samples"))
  
  return(counts)
}

#' Read sample information file
read_sample_info <- function(sample_info_file, count_matrix) {
  message("Reading sample information:", sample_info_file)
  
  # Read the sample information
  sample_info <- read.csv(sample_info_file)
  
  # Check if Sample column exists
  if (!"Sample" %in% colnames(sample_info)) {
    stop("Sample information must contain a 'Sample' column")
  }
  
  # Check if Group column exists
  if (!"Group" %in% colnames(sample_info)) {
    stop("Sample information must contain a 'Group' column")
  }
  
  # Set Sample as rownames
  rownames(sample_info) <- sample_info$Sample
  
  # Convert Group to factor
  sample_info$Group <- factor(sample_info$Group)
  
  # Check that sample names match between count matrix and sample info
  missing_samples <- setdiff(colnames(count_matrix), rownames(sample_info))
  if (length(missing_samples) > 0) {
    warning("Some samples in count matrix are not in sample information: ", 
            paste(missing_samples, collapse = ", "))
  }
  
  # Keep only samples that exist in the count matrix
  sample_info <- sample_info[intersect(rownames(sample_info), colnames(count_matrix)), ]
  
  # Check that we have sample information
  if (nrow(sample_info) == 0) {
    stop("No matching samples found between count matrix and sample information")
  }
  
  message(paste("Loaded sample information for", nrow(sample_info), "samples"))
  
  return(sample_info)
}

#' Perform differential binding analysis
run_differential_analysis <- function(counts, sample_info, output_dir, num_cores = 6) {
  message("Performing differential binding analysis...")
  
  # Ensure count matrix has samples in the same order as sample_info
  counts <- counts[, rownames(sample_info)]
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = sample_info, 
                                design = ~ Group)
  
  # Filter out bins with low counts
  keep <- rowSums(counts(dds)) > 1
  dds <- dds[keep, ]
  message(paste("Retained", sum(keep), "bins after filtering"))
  
  # Register parallel backend
  register(MulticoreParam(workers = num_cores))
  message(paste("Using", num_cores, "CPU cores for analysis"))
  
  # Run DESeq2 with parallelization
  message("Running DESeq2...")
  dds <- DESeq(dds, parallel = TRUE)
  
  # Define all possible pairwise comparisons
  groups <- levels(sample_info$Group)
  all_comparisons <- combn(groups, 2, simplify = FALSE)
  comparisons <- list()
  
  # Create comparison specifications for DESeq2
  for (pair in all_comparisons) {
    comparisons[[length(comparisons) + 1]] <- c("Group", pair[1], pair[2])
    comparisons[[length(comparisons) + 1]] <- c("Group", pair[2], pair[1])
  }
  
  # Print all available comparisons
  message("Available comparisons:")
  for (i in seq_along(comparisons)) {
    comp_name <- paste0(comparisons[[i]][2], "_vs_", comparisons[[i]][3])
    message(paste0("  ", i, ": ", comp_name))
  }
  
  return(list(dds = dds, comparisons = comparisons))
}

#' Process a specific contrast
process_contrast <- function(dds, comparison, output_dir, txdb, orgDb) {
  comp_name <- paste0(comparison[2], "_vs_", comparison[3])
  message("Processing comparison: ", comp_name)
  
  # Get results for this comparison
  res <- results(dds, contrast = comparison)
  
  # Order by adjusted p-value
  resOrdered <- res[order(res$padj), ]
  
  # Add bin coordinates
  bin_coords <- t(sapply(rownames(resOrdered), function(x) {
    parts <- unlist(strsplit(x, "[:-]"))
    if (length(parts) >= 3) {
      return(c(parts[1], parts[2], parts[3]))
    } else {
      return(c(NA, NA, NA))
    }
  }))
  
  resOrdered$Chromosome <- bin_coords[, 1]
  resOrdered$Start <- as.numeric(bin_coords[, 2])
  resOrdered$End <- as.numeric(bin_coords[, 3])
  
  # Convert to data frame
  resDF <- as.data.frame(resOrdered)
  
  # Save results
  result_file <- file.path(output_dir, paste0("differential_bins_", comp_name, ".csv"))
  write.csv(resDF, file = result_file, row.names = TRUE)
  
  # Create MA plot
  pdf(file.path(output_dir, paste0("MA_plot_", comp_name, ".pdf")))
  plotMA(res, main = paste0("MA Plot: ", comp_name))
  dev.off()
  
  # Filter significant bins
  message("Identifying significant differential bins...")
  sig_bins <- resDF[!is.na(resDF$padj) & resDF$padj < 0.05, ]
  
  message(paste("Found", nrow(sig_bins), "significant bins (padj < 0.05)"))
  
  # Save significant bins
  sig_file <- file.path(output_dir, paste0("significant_bins_", comp_name, ".csv"))
  write.csv(sig_bins, file = sig_file, row.names = TRUE)
  
  if (nrow(sig_bins) > 0) {
    # Create GRanges object for annotation
    bins_gr <- GRanges(
      seqnames = sig_bins$Chromosome,
      ranges = IRanges(start = sig_bins$Start, end = sig_bins$End),
      strand = "*",
      log2FoldChange = sig_bins$log2FoldChange,
      padj = sig_bins$padj
    )
    
    # Annotate bins
    message("Annotating significant bins...")
    binAnno <- annotatePeak(bins_gr, tssRegion = c(-3000, 3000),
                           TxDb = txdb, annoDb = "org.Mm.eg.db")
    
    # Save annotation summary plots
    pdf(file.path(output_dir, paste0("annotation_summary_", comp_name, ".pdf")))
    plotAnnoPie(binAnno)
    plotDistToTSS(binAnno)
    dev.off()
    
    # Convert to data frame
    binAnnoDF <- as.data.frame(binAnno)
    
    # Create a Bin_ID column for merging
    binAnnoDF$Bin_ID <- paste0(binAnnoDF$seqnames, ":", binAnnoDF$start, "-", binAnnoDF$end)
    
    # Add Bin_ID to results
    sig_bins$Bin_ID <- rownames(sig_bins)
    
    # Merge annotation with results
    res_annotated <- merge(sig_bins, binAnnoDF, by = "Bin_ID", all.x = TRUE)
    
    # Save annotated results
    anno_file <- file.path(output_dir, paste0("annotated_bins_", comp_name, ".csv"))
    write.csv(res_annotated, file = anno_file, row.names = FALSE)
    
    # Create volcano plot
    message("Creating volcano plot...")
    
    # Prepare data for volcano plot
    sig_threshold <- 0.05
    fc_threshold <- 1.0
    
    # Add significance categories
    res_annotated$Significance <- "Not Significant"
    res_annotated$Significance[res_annotated$padj < sig_threshold & res_annotated$log2FoldChange > fc_threshold] <- "Upregulated"
    res_annotated$Significance[res_annotated$padj < sig_threshold & res_annotated$log2FoldChange < -fc_threshold] <- "Downregulated"
    
    # Select top genes for labeling
    res_annotated$label <- NA
    top_up <- res_annotated[res_annotated$Significance == "Upregulated", ]
    if (nrow(top_up) > 0) {
      top_up <- top_up[order(top_up$padj)[1:min(15, nrow(top_up))], ]
      res_annotated$label[match(top_up$Bin_ID, res_annotated$Bin_ID)] <- top_up$SYMBOL
    }
    
    top_down <- res_annotated[res_annotated$Significance == "Downregulated", ]
    if (nrow(top_down) > 0) {
      top_down <- top_down[order(top_down$padj)[1:min(15, nrow(top_down))], ]
      res_annotated$label[match(top_down$Bin_ID, res_annotated$Bin_ID)] <- top_down$SYMBOL
    }
    
    # Calculate negative log10 p-values
    res_annotated$negLogPadj <- -log10(pmax(res_annotated$padj, 1e-300))
    
    # Create volcano plot
    volcano_plot <- ggplot(res_annotated, aes(x = log2FoldChange, y = negLogPadj, color = Significance)) +
      geom_point(size = 1, alpha = 0.7) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray70")) +
      geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "gray40") +
      geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "gray40") +
      geom_text_repel(aes(label = label), max.overlaps = 20, size = 3, box.padding = 0.5, force = 10) +
      labs(title = paste("Differential Binding:", comp_name),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted p-value") +
      theme_minimal() +
      theme(legend.position = "right")
    
    # Save volcano plot
    ggsave(file.path(output_dir, paste0("volcano_plot_", comp_name, ".pdf")), 
           plot = volcano_plot, width = 10, height = 8)
    
    # Perform functional enrichment analysis if there are enough genes
    if (sum(!is.na(res_annotated$ENTREZID)) >= 10) {
      message("Performing functional enrichment analysis...")
      
      # GO analysis
      tryCatch({
        ego <- enrichGO(gene = na.omit(unique(res_annotated$ENTREZID)),
                        OrgDb = orgDb,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
        
        if (nrow(ego) > 0) {
          # Save GO results
          go_file <- file.path(output_dir, paste0("GO_enrichment_", comp_name, ".csv"))
          write.csv(as.data.frame(ego), file = go_file, row.names = FALSE)
          
          # Create GO plots
          pdf(file.path(output_dir, paste0("GO_dotplot_", comp_name, ".pdf")), width = 10, height = 8)
          print(dotplot(ego, showCategory = 15, title = paste0("GO Enrichment: ", comp_name)))
          dev.off()
          
          # Create GO bar plot
          pdf(file.path(output_dir, paste0("GO_barplot_", comp_name, ".pdf")), width = 10, height = 8)
          print(barplot(ego, showCategory = 15, title = paste0("GO Enrichment: ", comp_name)))
          dev.off()
        } else {
          message("No significant GO terms found")
        }
      }, error = function(e) {
        message("Error in GO enrichment analysis: ", e$message)
      })
      
      # KEGG analysis
      tryCatch({
        ekegg <- enrichKEGG(gene = na.omit(unique(res_annotated$ENTREZID)),
                           organism = "mmu",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)
        
        if (nrow(ekegg) > 0) {
          # Save KEGG results
          kegg_file <- file.path(output_dir, paste0("KEGG_enrichment_", comp_name, ".csv"))
          write.csv(as.data.frame(ekegg), file = kegg_file, row.names = FALSE)
          
          # Create KEGG plots
          pdf(file.path(output_dir, paste0("KEGG_dotplot_", comp_name, ".pdf")), width = 10, height = 8)
          print(dotplot(ekegg, showCategory = 15, title = paste0("KEGG Enrichment: ", comp_name)))
          dev.off()
        } else {
          message("No significant KEGG pathways found")
        }
      }, error = function(e) {
        message("Error in KEGG enrichment analysis: ", e$message)
      })
    } else {
      message("Too few genes with annotations for enrichment analysis")
    }
  } else {
    message("No significant bins to annotate")
  }
  
  return(res)
}

#' Generate sample correlation and PCA plots
create_sample_plots <- function(dds, output_dir) {
  message("Creating sample visualization plots...")
  
  # Perform variance-stabilizing transformation
  vsd <- vst(dds, blind = FALSE)
  
  # Principal Component Analysis (PCA)
  pcaData <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # Create PCA plot
  pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = Group, label = name)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 10) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of Split-DamID Samples") +
    theme_minimal()
  
  # Save PCA plot
  ggsave(file.path(output_dir, "PCA_plot.pdf"), plot = pca_plot, width = 8, height = 6)
  
  # Sample distance heatmap
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste0(vsd$Group, "-", rownames(colData(vsd)))
  colnames(sampleDistMatrix) <- paste0(vsd$Group, "-", rownames(colData(vsd)))
  
  # Create heatmap of sample distances
  pdf(file.path(output_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
           main = "Sample Distance Matrix")
  dev.off()
}

#' Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check command line arguments
  if (length(args) < 3) {
    stop("Usage: Rscript splitdamid_diff_analysis.R <count_matrix> <sample_info> <output_dir> [contrast_number]")
  }
  
  count_file <- args[1]
  sample_info_file <- args[2]
  output_dir <- args[3]
  
  # Optional contrast number argument
  contrast_num <- NULL
  if (length(args) >= 4) {
    contrast_num <- as.integer(args[4])
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Determine number of CPU cores for parallel processing
  num_cores <- min(6, parallel::detectCores() - 1)
  if (num_cores < 1) num_cores <- 1
  
  # Set up log file
  log_file <- file.path(output_dir, "splitdamid_diff_analysis.log")
  con <- file(log_file, "w")
  sink(con, type = "output")
  sink(con, type = "message")
  
  message("Split-DamID Differential Analysis and Peak Annotation")
  message("=====================================================")
  message("Count matrix: ", count_file)
  message("Sample info: ", sample_info_file)
  message("Output directory: ", output_dir)
  message("CPU cores: ", num_cores)
  
  # Read input files
  counts <- read_count_matrix(count_file)
  sample_info <- read_sample_info(sample_info_file, counts)
  
  # Subset counts to match sample info
  counts <- counts[, rownames(sample_info)]
  
  # Run differential analysis
  result <- run_differential_analysis(counts, sample_info, output_dir, num_cores)
  dds <- result$dds
  comparisons <- result$comparisons
  
  # Create sample correlation and PCA plots
  create_sample_plots(dds, output_dir)
  
  # Load TxDb and OrgDb for annotation
  message("Loading genome databases for annotation...")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  orgDb <- org.Mm.eg.db
  
  # Process specific contrast or all contrasts
  if (!is.null(contrast_num) && contrast_num >= 1 && contrast_num <= length(comparisons)) {
    message("Processing selected contrast #", contrast_num)
    process_contrast(dds, comparisons[[contrast_num]], output_dir, txdb, orgDb)
  } else {
    message("Processing all contrasts")
    for (i in seq_along(comparisons)) {
      process_contrast(dds, comparisons[[i]], output_dir, txdb, orgDb)
    }
  }
  
  # Save normalized counts for future analysis
  message("Saving normalized counts...")
  normalized_counts <- counts(dds, normalized = TRUE)
  normalized_file <- file.path(output_dir, "normalized_counts.csv")
  write.csv(normalized_counts, file = normalized_file)
  
  # Save the R session info for reproducibility
  session_info <- sessionInfo()
  writeLines(capture.output(print(session_info)), file.path(output_dir, "sessionInfo.txt"))
  
  # Save the DESeq2 object for future analysis
  save_file <- file.path(output_dir, "splitdamid_deseq_analysis.RData")
  save(dds, comparisons, counts, sample_info, file = save_file)
  
  message("Analysis complete! Results saved to: ", output_dir)
  
  # Reset sink
  sink(type = "output")
  sink(type = "message")
  close(con)
}

# Run the main function
tryCatch({
  main()
}, error = function(e) {
  message("Error in analysis: ", e$message)
  if (sink.number() > 0) {
    sink(type = "output")
    sink(type = "message")
  }
  quit(status = 1)
})
