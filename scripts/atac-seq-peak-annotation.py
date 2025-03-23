#!/usr/bin/env python3
"""
ATAC-seq Peak Annotation Script

This script takes peak BED files from ATAC-seq analysis and annotates
them with nearby genes and genomic features.

Usage:
    python annotate_peaks.py --input peaks.bed --genome mm10 --output annotated_peaks.tsv

Author: [Your Name]
Date: March 23, 2025
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pyranges import PyRanges
import genomepy


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Annotate ATAC-seq peaks with genomic features')
    
    parser.add_argument('--input', '-i', required=True,
                        help='Input peak file (BED format)')
    
    parser.add_argument('--genome', '-g', required=True,
                        help='Genome assembly (e.g., mm10, hg38)')
    
    parser.add_argument('--output', '-o', required=True,
                        help='Output file for annotated peaks')
    
    parser.add_argument('--promoter-window', type=int, default=2000,
                        help='Window size around TSS to define promoters (default: 2000)')
    
    parser.add_argument('--enhancer-window', type=int, default=50000,
                        help='Maximum distance to assign distal peaks to genes (default: 50000)')
    
    return parser.parse_args()


def load_peaks(peak_file):
    """Load peak file into PyRanges object."""
    try:
        # Check file format based on extension
        if peak_file.endswith('.narrowPeak') or peak_file.endswith('.broadPeak'):
            # MACS2/MACS3 format
            columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 
                       'Strand', 'SignalValue', 'pValue', 'qValue']
            if peak_file.endswith('.narrowPeak'):
                columns.append('Peak')
                
            peaks_df = pd.read_csv(peak_file, sep='\t', header=None, names=columns)
            
        else:
            # Assume standard BED format
            peaks_df = pd.read_csv(peak_file, sep='\t', header=None)
            
            if peaks_df.shape[1] >= 3:
                # Minimum BED format requires at least 3 columns
                column_names = ['Chromosome', 'Start', 'End']
                
                # Add additional column names depending on the file's column count
                if peaks_df.shape[1] >= 4:
                    column_names.append('Name')
                if peaks_df.shape[1] >= 5:
                    column_names.append('Score')
                if peaks_df.shape[1] >= 6:
                    column_names.append('Strand')
                
                # Rename columns and keep any additional unnamed columns
                for i in range(len(column_names), peaks_df.shape[1]):
                    column_names.append(f'Column{i+1}')
                
                peaks_df.columns = column_names
            else:
                raise ValueError(f"BED file must have at least 3 columns. Found: {peaks_df.shape[1]}")
        
        # Create PyRanges object
        peaks = PyRanges(peaks_df)
        print(f"Loaded {len(peaks)} peaks from {peak_file}")
        return peaks
    
    except Exception as e:
        print(f"Error loading peak file {peak_file}: {e}")
        sys.exit(1)


def load_gene_annotations(genome):
    """Load gene annotations for the specified genome."""
    try:
        # Check if genome is already downloaded, otherwise download it
        if not genomepy.is_genome_downloaded(genome):
            print(f"Downloading {genome} genome and annotations...")
            genomepy.install_genome(genome, annotation=True)
        
        # Get annotation file path
        annotation_file = genomepy.Genome(genome).annotation
        
        if annotation_file is None or not os.path.exists(annotation_file):
            raise ValueError(f"Could not find annotation file for {genome}")
        
        # Load gene annotations
        genes = PyRanges.read_gtf(annotation_file)
        
        # Filter for genes only
        genes_df = genes.df.query("Feature == 'gene'").copy()
        
        print(f"Loaded {len(genes_df)} gene annotations from {genome}")
        return PyRanges(genes_df)
    
    except Exception as e:
        print(f"Error loading gene annotations for {genome}: {e}")
        sys.exit(1)


def define_promoters(genes, window_size=2000):
    """Define promoter regions around transcription start sites."""
    promoters = genes.copy()
    
    # For forward strand genes, the promoter is [start - window, start + window]
    forward_mask = promoters.df['Strand'] == '+'
    promoters.df.loc[forward_mask, 'Start'] = promoters.df.loc[forward_mask, 'Start'] - window_size
    promoters.df.loc[forward_mask, 'End'] = promoters.df.loc[forward_mask, 'Start'] + 2 * window_size
    
    # For reverse strand genes, the promoter is [end - window, end + window]
    reverse_mask = promoters.df['Strand'] == '-'
    promoters.df.loc[reverse_mask, 'Start'] = promoters.df.loc[reverse_mask, 'End'] - window_size
    promoters.df.loc[reverse_mask, 'End'] = promoters.df.loc[reverse_mask, 'End'] + window_size
    
    # Ensure start positions are not negative
    promoters.df.loc[promoters.df['Start'] < 0, 'Start'] = 0
    
    return promoters


def annotate_peaks_with_features(peaks, genes, promoters, enhancer_window=50000):
    """Annotate peaks with genomic features (promoter, gene body, etc.)."""
    # Initialize annotation columns
    peaks_df = peaks.df.copy()
    peaks_df['Feature'] = 'Intergenic'
    peaks_df['Nearest_Gene'] = ''
    peaks_df['Distance_to_TSS'] = np.nan
    peaks_df['Gene_ID'] = ''
    
    # Find peaks that overlap with promoters
    promoter_overlaps = peaks.join(promoters, how="left")
    
    if len(promoter_overlaps) > 0:
        # For overlapping peaks, mark them as promoters
        for idx, row in promoter_overlaps.df.iterrows():
            peak_idx = peaks_df.index[
                (peaks_df['Chromosome'] == row['Chromosome']) & 
                (peaks_df['Start'] == row['Start']) & 
                (peaks_df['End'] == row['End'])
            ].tolist()
            
            if peak_idx:
                peaks_df.loc[peak_idx, 'Feature'] = 'Promoter'
                peaks_df.loc[peak_idx, 'Nearest_Gene'] = row.get('gene_name', row.get('gene_id', ''))
                peaks_df.loc[peak_idx, 'Gene_ID'] = row.get('gene_id', '')
                
                # Calculate distance to TSS
                if row['Strand'] == '+':
                    tss = row['Start']
                else:
                    tss = row['End']
                
                peak_center = (row['Start'] + row['End']) / 2
                peaks_df.loc[peak_idx, 'Distance_to_TSS'] = peak_center - tss
    
    # Find peaks that overlap with gene bodies (but not promoters)
    gene_body_overlaps = peaks.join(genes, how="left")
    
    if len(gene_body_overlaps) > 0:
        for idx, row in gene_body_overlaps.df.iterrows():
            peak_idx = peaks_df.index[
                (peaks_df['Chromosome'] == row['Chromosome']) & 
                (peaks_df['Start'] == row['Start']) & 
                (peaks_df['End'] == row['End']) &
                (peaks_df['Feature'] == 'Intergenic')  # Only update if not already annotated
            ].tolist()
            
            if peak_idx:
                peaks_df.loc[peak_idx, 'Feature'] = 'Gene_Body'
                peaks_df.loc[peak_idx, 'Nearest_Gene'] = row.get('gene_name', row.get('gene_id', ''))
                peaks_df.loc[peak_idx, 'Gene_ID'] = row.get('gene_id', '')
                
                # Calculate distance to TSS
                if row['Strand'] == '+':
                    tss = genes.df.loc[genes.df['gene_id'] == row['gene_id'], 'Start'].values[0]
                else:
                    tss = genes.df.loc[genes.df['gene_id'] == row['gene_id'], 'End'].values[0]
                
                peak_center = (row['Start'] + row['End']) / 2
                peaks_df.loc[peak_idx, 'Distance_to_TSS'] = peak_center - tss
    
    # For intergenic peaks, find the nearest gene within the enhancer window
    intergenic_peaks = peaks_df[peaks_df['Feature'] == 'Intergenic'].copy()
    
    if len(intergenic_peaks) > 0:
        # For each intergenic peak, find the nearest gene
        for idx, peak in intergenic_peaks.iterrows():
            peak_center = (peak['Start'] + peak['End']) / 2
            min_distance = float('inf')
            nearest_gene = None
            nearest_gene_id = None
            nearest_tss = None
            
            for _, gene in genes.df.iterrows():
                if gene['Chromosome'] != peak['Chromosome']:
                    continue
                
                # Calculate TSS position
                if gene['Strand'] == '+':
                    tss = gene['Start']
                else:
                    tss = gene['End']
                
                # Calculate distance to TSS
                distance = abs(peak_center - tss)
                
                if distance < min_distance and distance <= enhancer_window:
                    min_distance = distance
                    nearest_gene = gene.get('gene_name', gene.get('gene_id', ''))
                    nearest_gene_id = gene.get('gene_id', '')
                    nearest_tss = tss
            
            if nearest_gene is not None:
                peaks_df.loc[idx, 'Feature'] = 'Enhancer'
                peaks_df.loc[idx, 'Nearest_Gene'] = nearest_gene
                peaks_df.loc[idx, 'Gene_ID'] = nearest_gene_id
                peaks_df.loc[idx, 'Distance_to_TSS'] = peak_center - nearest_tss
    
    return peaks_df


def add_genomic_context(peaks_df, genome):
    """Add additional genomic context like CpG islands, conservation, etc."""
    # This is a placeholder for additional annotation
    # Depending on available data for the genome, you could add:
    # - CpG islands
    # - Conservation scores
    # - Repeat elements
    # - Chromatin state from public datasets
    
    print("Note: Additional genomic context annotation can be implemented based on available data")
    
    return peaks_df


def main():
    """Main function to run the annotation pipeline."""
    # Parse command-line arguments
    args = parse_arguments()
    
    # Load peak regions
    peaks = load_peaks(args.input)
    
    # Load gene annotations
    genes = load_gene_annotations(args.genome)
    
    # Define promoter regions
    promoters = define_promoters(genes, window_size=args.promoter_window)
    
    # Annotate peaks
    peaks_df = annotate_peaks_with_features(
        peaks, genes, promoters, enhancer_window=args.enhancer_window
    )
    
    # Add additional genomic context if available
    peaks_df = add_genomic_context(peaks_df, args.genome)
    
    # Save annotated peaks
    print(f"Saving annotated peaks to {args.output}")
    peaks_df.to_csv(args.output, sep='\t', index=False)
    
    # Print summary
    feature_counts = peaks_df['Feature'].value_counts()
    print("\nPeak annotation summary:")
    for feature, count in feature_counts.items():
        print(f"  {feature}: {count} peaks ({count/len(peaks_df)*100:.1f}%)")
    
    print("\nAnnotation complete!")


if __name__ == "__main__":
    main()