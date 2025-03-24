#!/usr/bin/env python3

"""
Prepare Split-DamID Binned Count Matrix for Differential Analysis

This script processes bin count files from Split-DamID workflow and prepares
a count matrix suitable for differential binding analysis in R. It builds on
the workflow from splitdamid-workflow.sh which already creates binned data
but does not consolidate it into a single matrix with proper metadata.

Usage:
    python prepare_splitdamid_bins.py --input splitdamid_results --output count_matrix.csv [--sample-info sample_info.csv]

Author: Ching-Fang Chang
Date: March 23, 2025
"""

import os
import sys
import argparse
import glob
import pandas as pd
import re

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Prepare Split-DamID binned count matrix for differential analysis')
    
    parser.add_argument('--input', '-i', required=True,
                        help='Input directory with Split-DamID results (output from splitdamid-workflow.sh)')
    
    parser.add_argument('--output', '-o', required=True,
                        help='Output file for count matrix (CSV format)')
    
    parser.add_argument('--sample-info', '-s',
                        help='Sample information file (if not provided, will be created based on sample names)')
    
    parser.add_argument('--bin-dir', default='analysis',
                        help='Subdirectory within input dir containing bin counts (default: analysis)')
    
    parser.add_argument('--bin-pattern', default='*_bin_counts.bed',
                        help='Pattern to match bin count files (default: *_bin_counts.bed)')
    
    return parser.parse_args()

def find_bin_count_files(input_dir, bin_dir, bin_pattern):
    """Find all bin count files in the specified directory."""
    bin_path = os.path.join(input_dir, bin_dir)
    if not os.path.exists(bin_path):
        print(f"Error: Bin directory {bin_path} not found!")
        sys.exit(1)
    
    bin_files = glob.glob(os.path.join(bin_path, bin_pattern))
    if not bin_files:
        print(f"Error: No bin count files found in {bin_path} with pattern {bin_pattern}")
        sys.exit(1)
    
    print(f"Found {len(bin_files)} bin count files")
    return sorted(bin_files)

def infer_sample_groups(bin_files):
    """Infer sample groups from bin file names."""
    sample_groups = {}
    
    for file_path in bin_files:
        # Extract sample name from file path
        sample_name = os.path.basename(file_path).replace('_bin_counts.bed', '')
        
        # Determine group and replicate based on sample name patterns
        if re.search(r'Gli3-Hand2-Dox-\d+', sample_name):
            group = 'Gli3-Hand2-Dox'
            replicate = re.search(r'Gli3-Hand2-Dox-(\d+)', sample_name).group(1)
        elif re.search(r'Gli3-Hand2-\d+', sample_name):
            group = 'Gli3-Hand2'
            replicate = re.search(r'Gli3-Hand2-(\d+)', sample_name).group(1)
        elif re.search(r'DAM-Dox-\d+', sample_name):
            group = 'DAM-Dox'
            replicate = re.search(r'DAM-Dox-(\d+)', sample_name).group(1)
        elif re.search(r'DAM-\d+', sample_name):
            group = 'DAM'
            replicate = re.search(r'DAM-(\d+)', sample_name).group(1)
        else:
            # For other sample naming patterns
            group = 'Unknown'
            replicate = '1'
        
        sample_groups[sample_name] = {'Group': group, 'Replicate': replicate}
    
    return sample_groups

def create_sample_info(sample_groups, output_file):
    """Create a sample information CSV file."""
    # Create DataFrame from sample groups dictionary
    sample_info = pd.DataFrame.from_dict(sample_groups, orient='index')
    sample_info.index.name = 'Sample'
    sample_info.reset_index(inplace=True)
    
    # Save to CSV
    sample_info.to_csv(output_file, index=False)
    print(f"Created sample information file: {output_file}")
    
    return sample_info

def process_bin_counts(bin_files, sample_groups, output_file):
    """Process bin count files and create a count matrix."""
    print("Creating count matrix...")
    
    # Initialize with the first file to get bin coordinates
    first_file = bin_files[0]
    first_sample = os.path.basename(first_file).replace('_bin_counts.bed', '')
    
    # Read the first count file
    counts = pd.read_csv(first_file, sep='\t', header=None)
    counts.columns = ['Chromosome', 'Start', 'End', first_sample]
    
    # Create a 'Bin_ID' column by combining 'Chromosome', 'Start', and 'End'
    counts['Bin_ID'] = counts['Chromosome'] + ':' + counts['Start'].astype(str) + '-' + counts['End'].astype(str)
    
    # Set Bin_ID as index for easier merging
    count_matrix = counts[['Bin_ID', first_sample]].set_index('Bin_ID')
    
    # Add counts for each remaining sample
    for file in bin_files[1:]:
        sample_name = os.path.basename(file).replace('_bin_counts.bed', '')
        
        # Read only the count column (4th column)
        sample_counts = pd.read_csv(file, sep='\t', header=None, usecols=[3], names=[sample_name])
        
        # Create Bin_ID for this sample
        bin_df = pd.read_csv(file, sep='\t', header=None, usecols=[0, 1, 2], names=['Chromosome', 'Start', 'End'])
        bin_df['Bin_ID'] = bin_df['Chromosome'] + ':' + bin_df['Start'].astype(str) + '-' + bin_df['End'].astype(str)
        
        # Set index for joining
        sample_df = pd.DataFrame({sample_name: sample_counts[sample_name].values, 'Bin_ID': bin_df['Bin_ID']})
        sample_df.set_index('Bin_ID', inplace=True)
        
        # Add to count matrix
        count_matrix = count_matrix.join(sample_df)
    
    # Reset index to get Bin_ID as a column
    count_matrix = count_matrix.reset_index()
    
    # Add bin coordinates back
    bin_coords = counts[['Bin_ID', 'Chromosome', 'Start', 'End']]
    count_matrix = count_matrix.merge(bin_coords, on='Bin_ID')
    
    # Reorder columns to have bin information first
    cols = ['Bin_ID', 'Chromosome', 'Start', 'End'] + [col for col in count_matrix.columns if col not in ['Bin_ID', 'Chromosome', 'Start', 'End']]
    count_matrix = count_matrix[cols]
    
    # Save count matrix
    count_matrix.to_csv(output_file, index=False)
    print(f"Count matrix saved to {output_file}")
    
    # Print some stats
    print(f"  - Matrix dimensions: {count_matrix.shape[0]} bins × {count_matrix.shape[1]-4} samples")
    
    return count_matrix

def main():
    """Main function to run the processing pipeline."""
    # Parse command-line arguments
    args = parse_arguments()
    
    # Find bin count files
    bin_files = find_bin_count_files(args.input, args.bin_dir, args.bin_pattern)
    
    # Infer sample groups from filenames
    sample_groups = infer_sample_groups(bin_files)
    
    # Create sample info file if not provided
    if args.sample_info:
        sample_info_file = args.sample_info
        if not os.path.exists(sample_info_file):
            sample_info = create_sample_info(sample_groups, sample_info_file)
        else:
            print(f"Using existing sample information file: {sample_info_file}")
            sample_info = pd.read_csv(sample_info_file)
    else:
        # Create default sample info file in the same directory as output
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        sample_info_file = os.path.join(output_dir if output_dir else '.', 'sample_info.csv')
        sample_info = create_sample_info(sample_groups, sample_info_file)
    
    # Process bin counts and create matrix
    count_matrix = process_bin_counts(bin_files, sample_groups, args.output)
    
    print("\nProcessing complete!")
    print(f"Sample information: {sample_info_file}")
    print(f"Count matrix: {args.output}")
    print("\nNext steps:")
    print("1. Run your differential binding analysis with:")
    print(f"   Rscript splitdamid_diffbind_analysis.R {args.output} {sample_info_file}")

if __name__ == "__main__":
    main()
