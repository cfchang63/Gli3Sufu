#!/usr/bin/env python3
"""
Extend reads to the next GATC site for Split-DamID analysis.

This script extends read fragments in BEDPE format to the next GATC site in both directions,
which is essential for Split-DamID analysis.

Usage:
    python extend_reads_to_GATC_revised.py input.bedpe gatc_sites.bed output.bedpe

Author: Claude
Date: March 23, 2025
"""

import sys
import os
import csv
from collections import defaultdict
import bisect

def load_gatc_sites(gatc_file):
    """
    Load GATC sites from a BED file into a dictionary organized by chromosome.
    
    Args:
        gatc_file (str): Path to the GATC sites BED file
        
    Returns:
        dict: Dictionary with chromosomes as keys and sorted lists of GATC site positions as values
    """
    print(f"Loading GATC sites from {gatc_file}...")
    gatc_sites = defaultdict(list)
    
    with open(gatc_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            
            # Store the start position of each GATC site
            gatc_sites[chrom].append(start)
    
    # Sort the GATC sites for each chromosome
    for chrom in gatc_sites:
        gatc_sites[chrom].sort()
    
    # Print some statistics
    total_sites = sum(len(sites) for sites in gatc_sites.values())
    print(f"Loaded {total_sites} GATC sites across {len(gatc_sites)} chromosomes")
    
    return gatc_sites

def find_next_gatc(gatc_sites, chrom, position, direction):
    """
    Find the next GATC site in the specified direction.
    
    Args:
        gatc_sites (dict): Dictionary with chromosome as key and sorted list of GATC positions as value
        chrom (str): Chromosome name
        position (int): Current position
        direction (str): Direction to search, "upstream" or "downstream"
        
    Returns:
        int: Position of the next GATC site, or None if not found
    """
    if chrom not in gatc_sites:
        return None
    
    sites = gatc_sites[chrom]
    
    if direction == "downstream":
        # Find the first site greater than or equal to the position
        index = bisect.bisect_left(sites, position)
        if index < len(sites):
            return sites[index]
        return None
    
    elif direction == "upstream":
        # Find the first site less than or equal to the position
        index = bisect.bisect_right(sites, position) - 1
        if index >= 0:
            return sites[index]
        return None
    
    return None

def extend_reads(input_bedpe, gatc_sites, output_bedpe):
    """
    Extend read pairs in BEDPE format to the nearest GATC sites.
    
    Args:
        input_bedpe (str): Path to input BEDPE file
        gatc_sites (dict): Dictionary with GATC sites by chromosome
        output_bedpe (str): Path to output BEDPE file
    """
    print(f"Extending reads from {input_bedpe} to nearest GATC sites...")
    
    # Track statistics
    total_reads = 0
    extended_reads = 0
    skipped_reads = 0
    
    with open(input_bedpe, 'r') as infile, open(output_bedpe, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        for row in reader:
            total_reads += 1
            
            # Check if we have a valid BEDPE format with at least 6 columns
            if len(row) < 6:
                skipped_reads += 1
                continue
            
            chrom1, start1, end1, chrom2, start2, end2 = row[:6]
            
            # Skip if reads are on different chromosomes
            if chrom1 != chrom2:
                skipped_reads += 1
                continue
            
            # Convert positions to integers
            start1, end1 = int(start1), int(end1)
            start2, end2 = int(start2), int(end2)
            
            # Determine the overall fragment boundaries
            fragment_start = min(start1, start2)
            fragment_end = max(end1, end2)
            
            # Find nearest GATC sites
            nearest_upstream = find_next_gatc(gatc_sites, chrom1, fragment_start, "upstream")
            nearest_downstream = find_next_gatc(gatc_sites, chrom1, fragment_end, "downstream")
            
            # Skip if we can't find GATC sites
            if nearest_upstream is None or nearest_downstream is None:
                skipped_reads += 1
                continue
            
            # Create extended fragment
            extended_start = nearest_upstream
            extended_end = nearest_downstream + 4  # Add 4 to include the full GATC site
            
            # Write the extended fragment to output
            output_row = [chrom1, str(extended_start), str(extended_end)]
            
            # Append the rest of the columns if they exist
            if len(row) > 6:
                output_row.extend(row[6:])
            
            writer.writerow(output_row)
            extended_reads += 1
            
            # Print progress every 1 million reads
            if total_reads % 1000000 == 0:
                print(f"Processed {total_reads} reads...")
    
    print(f"Extension complete. Processed {total_reads} reads, extended {extended_reads} reads, skipped {skipped_reads} reads.")

def main():
    # Check command line arguments
    if len(sys.argv) != 4:
        print("Usage: python extend_reads_to_GATC_revised.py input.bedpe gatc_sites.bed output.bedpe")
        sys.exit(1)
    
    input_bedpe = sys.argv[1]
    gatc_file = sys.argv[2]
    output_bedpe = sys.argv[3]
    
    # Check if input files exist
    if not os.path.exists(input_bedpe):
        print(f"Error: Input BEDPE file {input_bedpe} not found!")
        sys.exit(1)
    
    if not os.path.exists(gatc_file):
        print(f"Error: GATC sites file {gatc_file} not found!")
        sys.exit(1)
    
    # Load GATC sites
    gatc_sites = load_gatc_sites(gatc_file)
    
    # Extend reads
    extend_reads(input_bedpe, gatc_sites, output_bedpe)
    
    print(f"Successfully extended reads to GATC sites. Output written to {output_bedpe}")

if __name__ == "__main__":
    main()