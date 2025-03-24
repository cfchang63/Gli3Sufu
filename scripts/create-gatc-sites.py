#!/usr/bin/env python3
"""
Create a BED file of GATC sites in a genome.

This script identifies all GATC restriction sites in a genome FASTA file
and outputs them as a BED file, which is required for Split-DamID analysis.

Usage:
    python create_gatc_sites.py genome.fa output.bed

Author: Ching-Fang Chang
Date: March 23, 2025
"""

import sys
import os
import re

def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a dictionary of sequences.
    
    Args:
        fasta_file (str): Path to the FASTA file
        
    Returns:
        dict: Dictionary with chromosome names as keys and sequences as values
    """
    sequences = {}
    current_chrom = None
    current_seq = []
    
    print(f"Reading genome from {fasta_file}...")
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Save the previous chromosome if it exists
                if current_chrom is not None:
                    sequences[current_chrom] = ''.join(current_seq)
                
                # Start a new chromosome
                header = line[1:].split()[0]  # Take the first word after '>'
                current_chrom = header
                current_seq = []
            else:
                current_seq.append(line.upper())
    
    # Save the last chromosome
    if current_chrom is not None:
        sequences[current_chrom] = ''.join(current_seq)
    
    print(f"Loaded {len(sequences)} chromosomes/contigs from the genome")
    return sequences

def find_gatc_sites(sequences, output_bed):
    """
    Find all GATC sites in the provided sequences and write them to a BED file.
    
    Args:
        sequences (dict): Dictionary with chromosome names as keys and sequences as values
        output_bed (str): Path to the output BED file
    """
    print(f"Finding GATC sites and writing to {output_bed}...")
    
    pattern = re.compile('GATC')
    total_sites = 0
    
    with open(output_bed, 'w') as f:
        for chrom, sequence in sequences.items():
            # Find all occurrences of GATC
            for match in pattern.finditer(sequence):
                start = match.start()
                end = match.end()
                f.write(f"{chrom}\t{start}\t{end}\tGATC\n")
                total_sites += 1
            
            print(f"Found {pattern.findall(sequence).count('GATC')} GATC sites in {chrom}")
    
    print(f"Total of {total_sites} GATC sites found and written to {output_bed}")

def main():
    # Check command line arguments
    if len(sys.argv) != 3:
        print("Usage: python create_gatc_sites.py genome.fa output.bed")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    output_bed = sys.argv[2]
    
    # Check if input file exists
    if not os.path.exists(fasta_file):
        print(f"Error: Genome FASTA file {fasta_file} not found!")
        sys.exit(1)
    
    # Parse the FASTA file
    sequences = parse_fasta(fasta_file)
    
    # Find GATC sites and write to BED file
    find_gatc_sites(sequences, output_bed)
    
    print(f"Successfully created GATC sites BED file at {output_bed}")

if __name__ == "__main__":
    main()
