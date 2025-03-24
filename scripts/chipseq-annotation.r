#!/usr/bin/env Rscript

#' ChIP-seq Peak Annotation
#' 
#' This script loads peak files from ChIP-seq experiments and performs 
#' comprehensive annotation to identify genomic features associated with peaks.
#' 
#' Usage:
#'   Rscript annotate_chipseq_peaks.R --input peak_file.bed --output annotated_peaks.txt [OPTIONS]
#'
#' Author: Ching-Fang Chang
#' Date: March 23, 2025
#' 
#' Usage:
#' Rscript chipseq-annotation.R --input peaks.bed --output annotated_peaks.txt --genome mm10 --tss-region -3000,3000 --go-analysis

# Load required libraries and handle installation if needed
suppressPackageStartupMessages({
  packages <- c("optparse", "ChIPseeker", "GenomicRanges", "TxDb.Mmusculus.UCSC.mm10.knownGene", 
                "org.Mm.eg.db", "biomaRt", "ggplot2", "dplyr", "stringr", "clusterProfiler")
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste0("Installing package: ", pkg))
      if (pkg %in% c("ChIPseeker", "GenomicRanges", "TxDb.Mmusculus.UCSC.mm10.knownGene", 
                     "org.Mm.eg.db", "clusterProfiler")) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager", repos = "http://cran.us.r-project.org")
        BiocManager::install(pkg, update = FALSE)
      } else {
        install.packages(pkg, repos = "http://cran.us.r-project.org")
      }
    }
    library(pkg, character.only = TRUE)
  }
})

# Parse command line arguments
option_list <- list(
  make_option(c("--input", "-i"), type = "character", default = NULL,
              help = "Input peak file in BED format"),
  make_option(c("--output", "-o"), type = "character", default = "annotated_peaks.txt",
              help = "Output file name for annotated peaks [default=%default]"),
  make_option(c("--genome", "-g"), type = "character", default = "mm10",
              help = "Genome assembly (mm10, hg38, etc.) [default=%default]"),
  make_option(c("--tss-region"), type = "character", default = "-3000,3000",
              help = "TSS region to define promoters (-3000,3000 means 3kb upstream and downstream) [default=%default]"),
  make_option(c("--flank-distance"), type = "integer", default = 5000,
              help = "Flank distance for gene assignment [default=%default]"),
  make_option(c("--plots-dir"), type = "character", default = "annotation_plots",
              help = "Directory for saving annotation plots [default=%default]"),
  make_option(c("--plot-format"), type = "character", default = "pdf",
              help = "Format for output plots (pdf, png, svg) [default=%default]"),
  make_option(c("--go-analysis"), type = "logical", default = TRUE, action = "store_true",
              help = "Perform GO enrichment analysis [default=%default]"),
  make_option(c("--skip-header"), type = "logical", default = FALSE, action = "store_true",
              help = "Skip header in BED file [default=%default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$input)) {
  stop("Input file is required. Use --input or -i option.")
}

# Create plots directory if needed
if (!dir.exists(opt$plots_dir)) {
  dir.create(opt$plots_dir, recursive = TRUE)
}

# Process TSS region parameter
tss_region <- as.numeric(strsplit(opt$tss_region, ",")[[1]])
if (length(tss_region) != 2) {
  stop("TSS region should be specified as two comma-separated numbers, e.g., -3000,3000")
}

# Load appropriate TxDb object based on genome
txdb <- NULL
org_db <- NULL

if (opt$genome == "mm10") {
  message("Loading mouse genome (mm10) annotation databases...")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org_db <- "org.Mm.eg.db"
} else if (opt$genome == "hg38") {
  message("Loading human genome (hg38) annotation databases...")
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", update = FALSE)
  }
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org_db <- "org.Hs.eg.db"
} else {
  stop(paste0("Unsupported genome: ", opt$genome, ". Currently supported: mm10, hg38"))
}

# Load BED file
message(paste0("Loading peak file: ", opt$input))
bed_file <- read.table(opt$input, 
                      header = FALSE, 
                      sep = "\t", 
                      skip = ifelse(opt$skip_header, 1, 0))

# Check if the file has the minimum required columns
if (ncol(bed_file) < 3) {
  stop("BED file must have at least 3 columns (chr, start, end)")
}

# Rename columns based on BED format
colnames_bed <- c("chr", "start", "end")
if (ncol(bed_file) >= 4) colnames_bed <- c(colnames_bed, "name")
if (ncol(bed_file) >= 5) colnames_bed <- c(colnames_bed, "score")
if (ncol(bed_file) >= 6) colnames_bed <- c(colnames_bed, "strand")
# Add generic names for any additional columns
if (ncol(bed_file) > 6) {
  for (i in 7:ncol(bed_file)) {
    colnames_bed <- c(colnames_bed, paste0("V", i))
  }
}
colnames(bed_file) <- colnames_bed[1:ncol(bed_file)]

# Convert BED to GRanges
message("Converting BED file to GRanges...")
if ("strand" %in% colnames(bed_file)) {
  peaks <- GRanges(
    seqnames = bed_file$chr,
    ranges = IRanges(start = bed_file$start, end = bed_file$end),
    strand = bed_file$strand
  )
} else {
  peaks <- GRanges(
    seqnames = bed_file$chr,
    ranges = IRanges(start = bed_file$start, end = bed_file$end),
    strand = "*"
  )
}

# Add name if available
if ("name" %in% colnames(bed_file)) {
  mcols(peaks)$name <- bed_file$name
}

# Add score if available
if ("score" %in% colnames(bed_file)) {
  mcols(peaks)$score <- bed_file$score
}

# Make sure sequence levels match the TxDb object
message("Ensuring sequence level compatibility...")
seqlevelsStyle(peaks) <- "UCSC"

# Filter out peaks on non-standard chromosomes
standard_chromosomes <- paste0("chr", c(1:19, "X", "Y", "M"))
peaks <- keepSeqlevels(peaks, value = intersect(seqlevels(peaks), standard_chromosomes), pruning.mode = "coarse")

# Annotate peaks
message("Annotating peaks...")
annotations <- annotatePeak(
  peaks,
  tssRegion = tss_region,
  TxDb = txdb,
  level = "transcript",
  assignGenomicAnnotation = TRUE,
  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                              "Downstream", "Intergenic"),
  annoDb = org_db,
  addFlankGeneInfo = TRUE,
  flankDistance = opt$flank_distance,
  sameStrand = FALSE,
  ignoreOverlap = FALSE,
  ignoreUpstream = FALSE,
  ignoreDownstream = FALSE,
  overlap = "TSS",
  verbose = TRUE
)

# Create summary plots
message("Generating annotation summary plots...")
# Annotation pie chart
pdf(file.path(opt$plots_dir, "peak_annotation_pie.pdf"), width = 10, height = 8)
plotAnnoPie(annotations)
dev.off()

# Annotation bar chart
pdf(file.path(opt$plots_dir, "peak_annotation_bar.pdf"), width = 10, height = 6)
plotAnnoBar(annotations)
dev.off()

# Distance to TSS
pdf(file.path(opt$plots_dir, "peak_distance_to_TSS.pdf"), width = 10, height = 6)
plotDistToTSS(annotations)
dev.off()

# Get annotation as data frame
annotations_df <- as.data.frame(annotations)

# Save full annotation data
message(paste0("Saving annotated peaks to: ", opt$output))
write.table(annotations_df, file = opt$output, sep = "\t", row.names = FALSE, quote = FALSE)

# Run GO enrichment analysis if requested
if (opt$go_analysis) {
  message("Performing GO enrichment analysis...")
  
  # Extract gene IDs
  gene_ids <- annotations_df$geneId
  gene_ids <- unique(gene_ids[!is.na(gene_ids)])
  
  if (length(gene_ids) > 0) {
    # GO Biological Process enrichment
    go_bp <- enrichGO(
      gene = gene_ids,
      OrgDb = org_db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    
    if (!is.null(go_bp) && nrow(go_bp) > 0) {
      # Save results
      write.table(as.data.frame(go_bp), 
                  file = paste0(tools::file_path_sans_ext(opt$output), "_GO_BP.txt"), 
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Plot top results
      pdf(file.path(opt$plots_dir, "GO_BP_dotplot.pdf"), width = 12, height = 8)
      print(dotplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment"))
      dev.off()
      
      pdf(file.path(opt$plots_dir, "GO_BP_enrichment.pdf"), width = 12, height = 10)
      print(enrichMap(go_bp, n = 30, vertex.label.cex = 0.7))
      dev.off()
    } else {
      message("No significant GO Biological Process terms found.")
    }
    
    # KEGG pathway enrichment
    kegg_pathway <- enrichKEGG(
      gene = gene_ids,
      organism = ifelse(opt$genome == "mm10", "mmu", "hsa"),
      keyType = "kegg",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
    
    if (!is.null(kegg_pathway) && nrow(kegg_pathway) > 0) {
      # Save results
      write.table(as.data.frame(kegg_pathway), 
                  file = paste0(tools::file_path_sans_ext(opt$output), "_KEGG.txt"), 
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Plot top results
      pdf(file.path(opt$plots_dir, "KEGG_pathway_dotplot.pdf"), width = 12, height = 8)
      print(dotplot(kegg_pathway, showCategory = 15, title = "KEGG Pathway Enrichment"))
      dev.off()
    } else {
      message("No significant KEGG pathways found.")
    }
  } else {
    message("No gene IDs found for enrichment analysis.")
  }
}

# Try to connect to biomaRt for additional gene information
message("Retrieving additional gene information from biomaRt...")
tryCatch({
  # Set up biomaRt connection
  ensembl <- useEnsembl(
    biomart = "genes", 
    dataset = ifelse(opt$genome == "mm10", "mmusculus_gene_ensembl", "hsapiens_gene_ensembl")
  )
  
  # Get unique gene IDs from annotations
  entrez_ids <- unique(annotations_df$geneId[!is.na(annotations_df$geneId)])
  
  if (length(entrez_ids) > 0) {
    # Get gene info from biomaRt
    gene_info <- getBM(
      attributes = c("entrezgene_id", "external_gene_name", "description", "chromosome_name", 
                    "start_position", "end_position", "strand", "gene_biotype"),
      filters = "entrezgene_id",
      values = entrez_ids,
      mart = ensembl
    )
    
    # Save gene info
    write.table(gene_info, 
                file = paste0(tools::file_path_sans_ext(opt$output), "_gene_info.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    message("Additional gene information saved.")
  }
}, error = function(e) {
  message("Warning: Could not retrieve additional gene information from biomaRt.")
  message(e$message)
})

message("Analysis complete!")

# Print summary information
cat("\nSummary of results:\n")
cat("------------------\n")
cat("Input file:", opt$input, "\n")
cat("Output file:", opt$output, "\n")
cat("Genome assembly:", opt$genome, "\n")
cat("Number of peaks analyzed:", length(peaks), "\n")
cat("Number of annotated peaks:", nrow(annotations_df), "\n")
cat("Number of unique genes:", length(unique(annotations_df$geneId[!is.na(annotations_df$geneId)])), "\n")
cat("Plots saved to:", opt$plots_dir, "\n\n")
