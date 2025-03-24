#!/usr/bin/env python3
"""
ATAC-seq/ChIP-seq/CUT&RUN Peak Enrichment Analysis

This script analyzes the enrichment of peaks in genomic features and gene ontology
terms. It works with the annotated peak files from HOMER's annotatePeaks.pl.

Usage:
    python analyze_peak_enrichment.py --input annotated_peaks.txt --output results_directory

Author: Ching-Fang Chang
Date: March 23, 2025
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze peak enrichment in genomic features')
    
    parser.add_argument('--input', '-i', required=True,
                        help='Input annotated peak file from HOMER')
    
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory for results')
    
    parser.add_argument('--genome', '-g', default='mm10',
                        help='Genome assembly (default: mm10)')
    
    parser.add_argument('--promoter-window', type=int, default=2000,
                        help='Window size around TSS to define promoters (default: 2000)')
    
    parser.add_argument('--background', '-b', 
                        help='Background peak file for enrichment analysis (optional)')
    
    parser.add_argument('--gene-list', '-l',
                        help='List of genes of interest for enrichment analysis (optional)')
    
    return parser.parse_args()

def load_annotated_peaks(file_path):
    """Load HOMER annotated peak file."""
    try:
        # HOMER's annotatePeaks.pl output has a variable number of columns
        # We'll read it and determine the format based on the header
        with open(file_path, 'r') as f:
            header = f.readline().strip().split('\t')
        
        # Read the file with pandas
        peaks_df = pd.read_csv(file_path, sep='\t', skiprows=0, header=0)
        
        # Check if this is an intersection file with multiple peak regions
        if len(header) > 20 and 'Distance to TSS' in header:
            print(f"Detected HOMER annotated peak file with {len(peaks_df)} peaks")
        else:
            print(f"File format not recognized as HOMER annotated peaks")
            print(f"Header: {header[:5]}...")
            sys.exit(1)
        
        return peaks_df
    
    except Exception as e:
        print(f"Error loading annotated peak file {file_path}: {e}")
        sys.exit(1)

def analyze_genomic_distribution(peaks_df, output_dir):
    """Analyze and visualize the genomic distribution of peaks."""
    print("Analyzing genomic distribution of peaks...")
    
    # Extract annotation categories
    annotation_col = 'Annotation' if 'Annotation' in peaks_df.columns else 'Feature'
    
    if annotation_col not in peaks_df.columns:
        print(f"Warning: Could not find annotation column. Available columns: {peaks_df.columns.tolist()}")
        return
    
    # Count occurrences of each gene
    gene_counts = peaks_df[gene_col].value_counts()
    
    # Save the top genes to a file
    top_genes = gene_counts.head(100)
    top_genes.to_csv(os.path.join(output_dir, 'top_genes.tsv'), sep='\t', header=True)
    
    # Plot top 20 genes
    plt.figure(figsize=(12, 8))
    top_20 = gene_counts.head(20)
    sns.barplot(x=top_20.values, y=top_20.index)
    plt.title('Top 20 Genes with Most Peaks')
    plt.xlabel('Number of Peaks')
    plt.ylabel('Gene Name')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'top_genes.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # If a gene list is provided, check for enrichment
    if gene_list is not None:
        try:
            with open(gene_list, 'r') as f:
                interest_genes = [line.strip() for line in f if line.strip()]
            
            # Check if genes of interest are present in peaks
            interest_genes_present = [gene for gene in interest_genes if gene in peaks_df[gene_col].values]
            
            # Save results
            with open(os.path.join(output_dir, 'gene_list_overlap.txt'), 'w') as f:
                f.write(f"Total genes in list: {len(interest_genes)}\n")
                f.write(f"Genes found in peaks: {len(interest_genes_present)}\n")
                f.write(f"Percentage overlap: {len(interest_genes_present)/len(interest_genes)*100:.2f}%\n\n")
                f.write("Genes found in peaks:\n")
                for gene in interest_genes_present:
                    f.write(f"{gene}\n")
            
            print(f"Gene list analysis complete. Found {len(interest_genes_present)} out of {len(interest_genes)} genes.")
        except Exception as e:
            print(f"Error processing gene list: {e}")
    
    # If background peaks are provided, perform enrichment analysis
    if background_df is not None:
        print("Performing enrichment analysis against background peaks...")
        
        # Get all genes in background
        background_genes = set(background_df[gene_col].dropna())
        target_genes = set(peaks_df[gene_col].dropna())
        
        # Get all unique genes
        all_genes = background_genes.union(target_genes)
        
        # Calculate enrichment for each gene
        enrichment_results = []
        
        for gene in all_genes:
            # Construct contingency table
            in_target = gene in target_genes
            in_background = gene in background_genes
            
            a = 1 if in_target else 0  # In target
            b = 1 if in_background else 0  # In background
            c = len(target_genes) - a  # Not in target
            d = len(background_genes) - b  # Not in background
            
            # Perform Fisher's exact test
            contingency_table = [[a, b], [c, d]]
            odds_ratio, p_value = fisher_exact(contingency_table)
            
            enrichment_results.append({
                'gene': gene,
                'in_target': in_target,
                'in_background': in_background,
                'odds_ratio': odds_ratio,
                'p_value': p_value
            })
        
        # Create DataFrame
        enrichment_df = pd.DataFrame(enrichment_results)
        
        # Apply multiple testing correction
        enrichment_df['adjusted_p_value'] = multipletests(
            enrichment_df['p_value'], method='fdr_bh'
        )[1]
        
        # Sort by p-value
        enrichment_df = enrichment_df.sort_values('p_value')
        
        # Save results
        enrichment_df.to_csv(os.path.join(output_dir, 'gene_enrichment.tsv'), sep='\t', index=False)
        
        # Plot top enriched genes
        plt.figure(figsize=(12, 8))
        top_enriched = enrichment_df.head(20)
        sns.barplot(x='-log10(p_value)', y='gene', data=top_enriched.assign(**{
            '-log10(p_value)': -np.log10(top_enriched['p_value'])
        }))
        plt.title('Top 20 Enriched Genes')
        plt.xlabel('-log10(p-value)')
        plt.ylabel('Gene Name')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'enriched_genes.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Enrichment analysis complete. Results saved to {output_dir}")

def create_bed_for_go_analysis(peaks_df, output_dir):
    """Create a BED file for GO analysis with GREAT."""
    print("Creating BED file for GO analysis...")
    
    # Extract required columns
    try:
        bed_df = peaks_df[['Chr', 'Start', 'End']].copy()
        
        # Write to BED file
        bed_file = os.path.join(output_dir, 'peaks_for_go.bed')
        bed_df.to_csv(bed_file, sep='\t', header=False, index=False)
        
        print(f"BED file for GO analysis created: {bed_file}")
        print("You can upload this file to GREAT (http://great.stanford.edu/public/html/) for GO enrichment analysis.")
    except Exception as e:
        print(f"Error creating BED file: {e}")

def main():
    """Main function to run the analysis pipeline."""
    # Parse command-line arguments
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Load annotated peaks
    peaks_df = load_annotated_peaks(args.input)
    
    # Load background peaks if provided
    background_df = None
    if args.background:
        print(f"Loading background peaks from {args.background}")
        background_df = load_annotated_peaks(args.background)
    
    # Analyze genomic distribution
    analyze_genomic_distribution(peaks_df, args.output)
    
    # Analyze distance to TSS
    analyze_distance_to_tss(peaks_df, args.output, window=args.promoter_window)
    
    # Analyze gene enrichment
    analyze_gene_enrichment(peaks_df, args.output, background_df, args.gene_list)
    
    # Create BED file for GO analysis
    create_bed_for_go_analysis(peaks_df, args.output)
    
    print(f"\nAnalysis complete! Results saved to {args.output}")

if __name__ == "__main__":
    main()urrences of each annotation category
    annotation_counts = peaks_df[annotation_col].value_counts()
    
    # Create a simplified category for visualization
    def simplify_category(annotation):
        if 'promoter' in annotation.lower() or 'tss' in annotation.lower():
            return 'Promoter'
        elif 'exon' in annotation.lower():
            return 'Exon'
        elif 'intron' in annotation.lower():
            return 'Intron'
        elif '3\'' in annotation.lower() or '3 utr' in annotation.lower():
            return '3\' UTR'
        elif '5\'' in annotation.lower() or '5 utr' in annotation.lower():
            return '5\' UTR'
        elif 'tts' in annotation.lower():
            return 'Transcription Termination Site'
        elif 'intergenic' in annotation.lower():
            return 'Intergenic'
        else:
            return 'Other'
    
    peaks_df['Simple_Category'] = peaks_df[annotation_col].apply(simplify_category)
    simple_counts = peaks_df['Simple_Category'].value_counts()
    
    # Save counts to file
    annotation_counts.to_csv(os.path.join(output_dir, 'annotation_distribution.tsv'), sep='\t', header=True)
    simple_counts.to_csv(os.path.join(output_dir, 'simplified_distribution.tsv'), sep='\t', header=True)
    
    # Create pie charts
    plt.figure(figsize=(12, 5))
    
    # Detailed pie chart
    plt.subplot(1, 2, 1)
    plt.pie(annotation_counts.values, labels=annotation_counts.index, autopct='%1.1f%%')
    plt.title('Detailed Genomic Distribution')
    
    # Simplified pie chart
    plt.subplot(1, 2, 2)
    plt.pie(simple_counts.values, labels=simple_counts.index, autopct='%1.1f%%')
    plt.title('Simplified Genomic Distribution')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'genomic_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Genomic distribution analysis complete. Results saved to {output_dir}")

def analyze_distance_to_tss(peaks_df, output_dir, window=2000):
    """Analyze and visualize peak distance to TSS."""
    print("Analyzing peak distance to TSS...")
    
    distance_col = 'Distance to TSS' if 'Distance to TSS' in peaks_df.columns else None
    
    if distance_col is None:
        # Try to find a column that contains distance information
        for col in peaks_df.columns:
            if 'distance' in col.lower() and 'tss' in col.lower():
                distance_col = col
                break
    
    if distance_col is None:
        print("Warning: Could not find distance to TSS column")
        return
    
    # Convert distance to numeric and handle non-numeric values
    peaks_df[distance_col] = pd.to_numeric(peaks_df[distance_col], errors='coerce')
    
    # Filter out NaN values
    distance_df = peaks_df[~peaks_df[distance_col].isna()]
    
    # Create histogram
    plt.figure(figsize=(10, 6))
    
    # Filter distances to a reasonable range for visualization
    distances = distance_df[distance_col]
    filtered_distances = distances[(distances >= -window) & (distances <= window)]
    
    # Create histogram
    sns.histplot(filtered_distances, bins=50, kde=True)
    plt.axvline(x=0, color='r', linestyle='--', alpha=0.7)
    plt.title(f'Distribution of Peak Distance to TSS (±{window} bp)')
    plt.xlabel('Distance to TSS (bp)')
    plt.ylabel('Number of Peaks')
    plt.grid(alpha=0.3)
    
    # Save plot
    plt.savefig(os.path.join(output_dir, 'distance_to_tss.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Calculate statistics
    stats = {
        'total_peaks': len(distance_df),
        'mean_distance': distances.mean(),
        'median_distance': distances.median(),
        'promoter_peaks': len(filtered_distances),
        'percent_promoter': len(filtered_distances) / len(distance_df) * 100 if len(distance_df) > 0 else 0
    }
    
    # Save statistics
    with open(os.path.join(output_dir, 'tss_distance_stats.txt'), 'w') as f:
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")
    
    print(f"TSS distance analysis complete. Results saved to {output_dir}")

def analyze_gene_enrichment(peaks_df, output_dir, background_df=None, gene_list=None):
    """Analyze and visualize gene enrichment in the peaks."""
    print("Analyzing gene enrichment...")
    
    # Extract gene names from annotated peaks
    gene_col = 'Nearest Promoter' if 'Nearest Promoter' in peaks_df.columns else 'Nearest_Gene'
    
    if gene_col not in peaks_df.columns:
        # Try to find a suitable column
        for col in peaks_df.columns:
            if 'gene' in col.lower() or 'promoter' in col.lower():
                gene_col = col
                break
    
    if gene_col not in peaks_df.columns:
        print(f"Warning: Could not find gene name column. Available columns: {peaks_df.columns.tolist()}")
        return
    
    # Count occ