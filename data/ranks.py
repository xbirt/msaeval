import os
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
import re

# This script takes a list of reference sequences and their corresponding ranks as output by grinder
# groups the main taxonomies and their corresponding percentages and prints them out.

# Paths to input files (adjust if needed)
RANKS_PATH = os.path.join(os.path.dirname(__file__), '../500k-ranks.txt')
FASTA_PATH = os.path.join(os.path.dirname(__file__), '../selected_species_500000_subset.fasta')

def parse_ranks_file(ranks_path):
    ref_to_percent = {}
    with open(ranks_path, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 3:
                continue
            _, ref, percent = parts
            ref_to_percent[ref] = float(percent)
    return ref_to_percent

def parse_fasta_headers(fasta_path):
    ref_to_taxonomy = {}
    for record in SeqIO.parse(fasta_path, 'fasta'):
        header = record.description
        if ' ' in header:
            ref, taxonomy = header.split(' ', 1)
            ref_to_taxonomy[ref] = taxonomy
    return ref_to_taxonomy

def extract_main_taxonomy(taxonomy):
    # Match up to the first word after the last semicolon
    # Example: '...;Escherichia coli K-12' -> '...;Escherichia coli'
    match = re.match(r'^(.*;\s*\S+\s+\S+)', taxonomy)
    if match:
        return match.group(1)
    # Fallback: try up to the first two words after last semicolon
    match = re.match(r'^(.*;\s*\S+)', taxonomy)
    if match:
        return match.group(1)
    return taxonomy

def main():
    ref_to_percent = parse_ranks_file(RANKS_PATH)
    ref_to_taxonomy = parse_fasta_headers(FASTA_PATH)

    taxonomy_to_percent = defaultdict(float)
    for ref, percent in ref_to_percent.items():
        taxonomy = ref_to_taxonomy.get(ref)
        if taxonomy:
            main_taxonomy = extract_main_taxonomy(taxonomy)
            taxonomy_to_percent[main_taxonomy] += percent
        else:
            print(f'Warning: Reference {ref} not found in FASTA headers.')

    # Sort by percentage descending
    sorted_taxa = sorted(taxonomy_to_percent.items(), key=itemgetter(1), reverse=True)
    with open('ranks.txt', 'w') as out_f:
        for taxonomy, percent in sorted_taxa:
            last_part = taxonomy.split(';')[-1].strip()
            if last_part == 'uncultured bacterium':
                continue
            line = f'{taxonomy},{last_part},{percent}'
            print(line)
            out_f.write(line + '\n')

if __name__ == '__main__':
    main()
