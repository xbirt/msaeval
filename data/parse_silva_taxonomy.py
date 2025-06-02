#! /usr/bin/env python3

import csv
import os
from Bio import SeqIO

"""
SILVA Taxonomy Parser

This script parses SILVA reference taxonomy FASTA files and extracts taxonomic 
information for each sequence accession. It handles flexible taxonomy structures
for different domains (Bacteria, Archaea, Eukaryota) and saves the results to a
comprehensive CSV file.

Features:
- Flexible parsing of varying taxonomy levels
- Support for prokaryotic and eukaryotic taxonomy structures  
- Comprehensive field extraction including domain, phylum, class, order, family, genus, species
- Special handling for eukaryotic supergroups and divisions
- Progress tracking during file processing
- CSV export with all taxonomy fields

Output: silva_accession_taxa.csv with columns for accession and all taxonomy levels
"""

def parse_silva_taxonomy_flexible(taxonomy_string):
    """
    Parse Silva taxonomy string with flexible handling of varying levels.
    
    Examples:
    'Bacteria;Proteobacteria;Gammaproteobacteria;Burkholderiales;Comamonadaceae;Mitsuaria;unidentified'
    'Eukaryota;SAR;Alveolata;Dinoflagellata;Dinophyceae;Peridiniphycidae;Gonyaulacales;Alexandrium;Alexandrium tamarense'
    """
    if not taxonomy_string:
        return {}
    
    # Split by semicolon and clean
    taxa_levels = [level.strip() for level in taxonomy_string.split(';') if level.strip()]
    
    if not taxa_levels:
        return {}
    
    # Create a flexible taxonomy structure
    taxonomy_dict = {
        'full_lineage': taxa_levels,
        'lineage_length': len(taxa_levels),
        'domain': taxa_levels[0] if len(taxa_levels) > 0 else None,
        'lowest_level': taxa_levels[-1] if taxa_levels else None
    }
    
    # Try to identify standard ranks based on common patterns and position
    if len(taxa_levels) >= 1:
        taxonomy_dict['domain'] = taxa_levels[0]
    
    # For bacteria/archaea, try to map common positions
    if taxa_levels[0] in ['Bacteria', 'Archaea']:
        if len(taxa_levels) >= 2:
            taxonomy_dict['phylum'] = taxa_levels[1]
        if len(taxa_levels) >= 3:
            taxonomy_dict['class'] = taxa_levels[2]
        if len(taxa_levels) >= 4:
            taxonomy_dict['order'] = taxa_levels[3]
        if len(taxa_levels) >= 5:
            taxonomy_dict['family'] = taxa_levels[4]
        if len(taxa_levels) >= 6:
            if len(taxa_levels) == 6:
                # For 6 levels, split the last level into genus and species
                last_level = taxa_levels[5]
                parts = last_level.split(' ', 1)
                taxonomy_dict['genus'] = parts[0]
                taxonomy_dict['species'] = last_level
            else:
                taxonomy_dict['genus'] = taxa_levels[5]
        if len(taxa_levels) >= 7:
            taxonomy_dict['genus'] = taxa_levels[5]
            taxonomy_dict['species'] = taxa_levels[6]
    
    # For eukaryotes, the structure can be more variable
    elif taxa_levels[0] == 'Eukaryota':
        # Try to identify genus (usually second to last) and species (last)
        if len(taxa_levels) >= 2:
            # Try to find genus-like names (capitalized, might be followed by species)
            for i in range(len(taxa_levels) - 1, 0, -1):
                level = taxa_levels[i]
                # Check if this looks like a genus (starts with capital, not obviously higher rank)
                if level and level[0].isupper() and ' ' not in level:
                    taxonomy_dict['genus'] = level
                    if i + 1 < len(taxa_levels):
                        taxonomy_dict['species'] = taxa_levels[i + 1]
                    break
        
        # Try to assign higher level taxa for eukaryotes
        if len(taxa_levels) >= 2:
            taxonomy_dict['supergroup'] = taxa_levels[1]  # SAR, Archaeplastida, etc.
        if len(taxa_levels) >= 3:
            taxonomy_dict['division'] = taxa_levels[2]
    
    return taxonomy_dict

def load_silva_reference(silva_file="SILVA_138_SSURef_tax_silva.fasta"):
    """Load and parse SILVA reference file."""
    silva_taxonomy = {}
    record_count = 0
    
    try:
        print(f"Loading SILVA reference file: {silva_file}")
        for record in SeqIO.parse(silva_file, "fasta"):
            record_count += 1
            header = record.description
            
            # Parse header: >accession taxonomy_string
            parts = header.split(' ', 1)
            if len(parts) == 2:
                accession = parts[0]
                taxonomy_string = parts[1]
                silva_taxonomy[accession] = parse_silva_taxonomy_flexible(taxonomy_string)
            
            # Show progress every 1000 records
            if record_count % 1000 == 0:
                print(f"  Processed {record_count} records...")
        
        print(f"Loaded {len(silva_taxonomy)} SILVA reference taxonomies (total records processed: {record_count})")
    except Exception as e:
        print(f"Warning: Could not load SILVA reference file {silva_file}: {e}")
    
    return silva_taxonomy

def save_silva_taxonomy_to_csv(silva_taxonomy, output_file="silva_accession_taxa.csv"):
    """Save SILVA taxonomy data to CSV file with all possible fields."""
    
    # Define all possible fields we might encounter
    fieldnames = [
        'accession',
        'full_lineage',
        'lineage_length',
        'domain',
        'lowest_level',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'supergroup',  # For eukaryotes
        'division'     # For eukaryotes
    ]
    
    print(f"Saving taxonomy data to {output_file}")
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        record_count = 0
        for accession, taxonomy_data in silva_taxonomy.items():
            record_count += 1
            
            # Create row with all fields, using empty strings for missing values
            row = {'accession': accession}
            
            # Convert full_lineage list to semicolon-separated string
            if 'full_lineage' in taxonomy_data:
                row['full_lineage'] = ';'.join(taxonomy_data['full_lineage'])
            
            # Add all other fields
            for field in fieldnames[1:]:  # Skip accession as we already added it
                if field == 'full_lineage':
                    continue  # Already handled above
                row[field] = taxonomy_data.get(field, '')
        
            writer.writerow(row)
            
            # Show progress every 1000 records
            if record_count % 1000 == 0:
                print(f"  Saved {record_count} records...")
    
    print(f"Successfully saved {len(silva_taxonomy)} records to {output_file}")

def main():
    """Main function to parse SILVA FASTA file and save to CSV."""
    
    # Look for SILVA FASTA file in current directory
    silva_files = [
        "SILVA_138_SSURef_tax_silva.fasta",
        "SILVA_138_SSURef_NR99_tax_silva.fasta",
        "silva_reference.fasta"
    ]
    
    silva_file = None
    for filename in silva_files:
        if os.path.exists(filename):
            silva_file = filename
            break
    
    if not silva_file:
        print("No SILVA FASTA file found. Please ensure one of the following files exists:")
        for filename in silva_files:
            print(f"  - {filename}")
        print("\nOr specify the file path when calling load_silva_reference()")
        return
    
    # Load and parse SILVA taxonomy
    silva_taxonomy = load_silva_reference(silva_file)
    
    if not silva_taxonomy:
        print("No taxonomy data loaded. Please check the FASTA file format.")
        return
    
    # Save to CSV
    save_silva_taxonomy_to_csv(silva_taxonomy, "silva_accession_taxa.csv")
    
    # Print some statistics
    print(f"\nSummary:")
    print(f"Total records: {len(silva_taxonomy)}")
    
    # Count domains
    domains = {}
    for taxonomy_data in silva_taxonomy.values():
        domain = taxonomy_data.get('domain', 'Unknown')
        domains[domain] = domains.get(domain, 0) + 1
    
    print("Domain distribution:")
    for domain, count in sorted(domains.items()):
        print(f"  {domain}: {count}")

if __name__ == "__main__":
    main()