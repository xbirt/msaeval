#!/usr/bin/env python3

import os
import re
import glob
import csv
from pathlib import Path
from Bio import SeqIO

def load_mothur_silva_mappings():
    """Load mothur to silva reference mappings from TSV files."""
    seed_mapping = {}
    nr99_mapping = {}
    
    # Load seed mapping
    try:
        with open('mothur2silva-seed.txt', 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) >= 2:
                    mothur_ref = row[0].strip()
                    silva_ref = row[1].strip()
                    if mothur_ref and silva_ref:
                        seed_mapping[mothur_ref] = silva_ref
        print(f"Loaded {len(seed_mapping)} seed mappings from mothur2silva-seed.txt")
    except Exception as e:
        print(f"Warning: Could not load mothur2silva-seed.txt: {e}")
    
    # Load nr99 mapping
    try:
        with open('mothur2silva-nr99.txt', 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) >= 2:
                    mothur_ref = row[0].strip()
                    silva_ref = row[1].strip()
                    if mothur_ref and silva_ref:
                        nr99_mapping[mothur_ref] = silva_ref
        print(f"Loaded {len(nr99_mapping)} nr99 mappings from mothur2silva-nr99.txt")
    except Exception as e:
        print(f"Warning: Could not load mothur2silva-nr99.txt: {e}")
    
    return seed_mapping, nr99_mapping

def extract_command_info(command_line):
    """Extract executable and parameters from command line."""
    if not command_line.startswith("Comanda: "):
        return None, None
    
    command = command_line[9:].strip()  # Remove "Comanda: " prefix
    
    # Split into executable and parameters
    parts = command.split(' ', 1)
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], parts[1]

def determine_label(executable, parameters):
    """Determine label based on executable and parameters."""
    # Check for kraken2-build first (to disregard)
    if 'kraken2-build' in executable:
        return None
    
    executable = executable.lower()
    parameters = parameters.lower()

    # Determine base label
    if 'mothur' in executable:
        label = 'mothur'
    elif 'sina' in executable:
        label = 'sina'
    elif 'minimap2' in executable:
        label = 'minimap2'
    elif 'kraken2' in executable:
        label = 'kraken2'
    elif 'nast' in executable or 'pynast' in executable:
        label = 'nast'
    else:
        return None
    
    # Check for lowset suffix
    if 'seed_' in parameters:
        label += '-seed'
    elif 'nr99_' in parameters:
        label += '-nr99'
    elif 'nr_' in parameters:
        label += '-nr'
    return label

def find_results_file(label, parameters, log_file_dir):
    """Find the results file based on command type and parameters."""
    results_file = None
    
    if label.startswith('mothur'):
        # Find input fasta file with pattern ([0-9.]+[km])\.fasta
        fasta_match = re.search(r'([0-9.]+[km])\.fasta', parameters)
        if fasta_match:
            fasta_name = fasta_match.group(1)
            results_file = f"{fasta_name}.align_report"
    
    elif label.startswith('sina'):
        # Find CSV output file with pattern --out ([^ ]+\.csv)
        csv_match = re.search(r'--out ([^ ]+\.csv)', parameters)
        if csv_match:
            results_file = csv_match.group(1)
    
    elif label.startswith('minimap2'):
        # Find PAF output file with pattern > ([^ ]+\.paf)
        paf_match = re.search(r'> ([^ ]+\.paf)', parameters)
        if paf_match:
            results_file = paf_match.group(1)
    
    elif label.startswith('kraken2'):
        # Find report file with pattern --output ([^ ]+)
        report_match = re.search(r'--output ([^ ]+)', parameters)
        if report_match:
            results_file = report_match.group(1)
    
    return results_file

def extract_input_fasta(parameters):
    """Extract the input fasta file base name from parameters."""
    fasta_match = re.search(r'([0-9.]+[km])\.fasta', parameters)
    if fasta_match:
        return fasta_match.group(1)  # Return base name without extension
    return None

def parse_results_file(file_path, label, col1_idx, col2_idx, seed_mapping=None, nr99_mapping=None):
    """Parse results file and extract ID->mapping pairs."""
    id_mapping = {}
    
    try:
        # Determine if it's CSV or TSV
        is_csv = label.startswith('sina')
        delimiter = ',' if is_csv else '\t'
        
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter=delimiter)
            
            for row in reader:
                if len(row) > max(col1_idx, col2_idx):
                    read_id = row[col1_idx].strip()
                    mapping = row[col2_idx].strip()

                    if not read_id.isdigit():
                        continue

                    # Skip if ID already exists (keep first occurrence only)
                    if read_id in id_mapping:
                        continue
                    
                    # If mapping has multiple options, keep only the first
                    if ' ' in mapping:
                        mapping = mapping.split()[0]

                    # If mapping has a colon, keep only the part before it
                    if ':' in mapping:
                        mapping = mapping.split(':')[0]
                    
                    # Apply mothur to silva mapping if applicable
                    if label.startswith('mothur') and mapping:
                        if 'seed' in label and seed_mapping and mapping in seed_mapping:
                            mapping = seed_mapping[mapping]
                        elif 'nr' in label and nr99_mapping and mapping in nr99_mapping:
                            mapping = nr99_mapping[mapping]
                    
                    if read_id and mapping:
                        id_mapping[read_id] = mapping
    
    except Exception as e:
        print(f"Error parsing results file {file_path}: {e}")
    
    return id_mapping

def analyze_fasta_identity(fasta_file, id_mapping):
    """Analyze FASTA file and calculate identity based on mappings."""
    matches = 0
    total = 0
    
    try:
        if not os.path.exists(fasta_file):
            print(f"Warning: FASTA file {fasta_file} not found")
            return 0, 0
        
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description
            
            # Extract ID (number after > followed by space)
            id_match = re.match(r'^(\d+)\s', header)
            if not id_match:
                continue
            
            read_id = id_match.group(1)
            total += 1
            
            # Look up mapping for this ID
            if read_id in id_mapping:
                mapping = id_mapping[read_id].lower()
                header_lower = header.lower()
                
                # Check if mapping is found in header
                if mapping in header_lower:
                    matches += 1
    
    except Exception as e:
        print(f"Error analyzing FASTA file {fasta_file}: {e}")
    
    return matches, total

def get_column_indices(label):
    """Get column indices for different command types."""
    if label.startswith('sina'):
        return 0, 3  # columns 1 and 4 (0-indexed)
    elif label.startswith('kraken2'):
        return 1, 4  # columns 2 and 5 (0-indexed)
    elif label.startswith('minimap2'):
        return 0, 5  # columns 1 and 6 (0-indexed)
    elif label.startswith('mothur'):
        return 0, 2  # columns 1 and 3 (0-indexed)
    else:
        return None, None

def process_benchmark_files(base_folder):
    """Process all benchmark files in the base folder recursively."""
    results = []
    base_path = Path(base_folder)
    
    # Load mothur to silva mappings
    seed_mapping, nr99_mapping = load_mothur_silva_mappings()
    
    # Find all *_detalii.log files recursively
    pattern = "**/benchmark_*_detalii.log"
    for log_file in base_path.glob(pattern):
        try:
            # Extract timestamp from filename
            filename = log_file.name
            timestamp_match = re.search(r'benchmark_(\d+)_detalii\.log', filename)
            if not timestamp_match:
                continue
            
            timestamp = timestamp_match.group(1)
            
            # Read first line of the file
            with open(log_file, 'r', encoding='utf-8') as f:
                first_line = f.readline().strip()
            
            if not first_line:
                continue
            
            # Extract command info
            executable, parameters = extract_command_info(first_line)
            if not executable:
                continue
            
            # Determine label
            label = determine_label(executable, parameters)
            if not label:  # Skip if no valid label or kraken2-build
                continue
            
            # Extract input fasta file name
            input_fasta = extract_input_fasta(parameters)
            
            # Find results file
            results_file = find_results_file(label, parameters, log_file.parent)
            if not results_file:
                print(f"Error: Could not determine results file for {log_file}")
                continue
            
            # Check if results file exists
            results_path = log_file.parent / results_file
            if not results_path.exists():
                print(f"Error: Results file {results_file} does not exist in {log_file.parent}")
                continue
            
            if not label.startswith('sina'):
                #continue
                pass

            if 'ncbi' in str(results_path):
                continue
            
            # Parse results file and get ID->mapping pairs
            col1_idx, col2_idx = get_column_indices(label)
            if col1_idx is None or col2_idx is None:
                print(f"Error: Unknown column indices for label {label}")
                continue
            
            id_mapping = parse_results_file(results_path, label, col1_idx, col2_idx, seed_mapping, nr99_mapping)
            
            # Analyze FASTA file for identity
            total_reads = 0
            identical_matches = 0
            match_percent = "0.00"
            
            if input_fasta:
                fasta_file = f"{input_fasta}.fasta"
                matches, total = analyze_fasta_identity(fasta_file, id_mapping)
                
                total_reads = total
                identical_matches = matches
                
                if total > 0:
                    percentage = (matches / total) * 100
                    match_percent = f"{percentage:.2f}"
            
            # Get relative path from base folder
            relative_path = log_file.parent.relative_to(base_path)
            
            results.append({
                'relative_path': str(relative_path),
                'timestamp': timestamp,
                'label': label,
                'input_fasta': input_fasta or 'N/A',
                'results_file': results_file,
                'total_reads': str(total_reads),
                'identical_matches': str(identical_matches),
                'match_percent': match_percent + '%'
            })
            
        except Exception as e:
            print(f"Error processing {log_file}: {e}")
    
    return results

def print_ascii_table(results):
    """Print results in a nicely formatted ASCII table."""
    if not results:
        print("No valid benchmark files found.")
        return
    
    # Calculate column widths
    headers = ['File Path', 'Timestamp', 'Command Label', 'Input Fasta', 'Results File', 'Total Reads', 'Identical Matches', 'Match Percent']
    col_widths = [len(header) for header in headers]
    
    for result in results:
        col_widths[0] = max(col_widths[0], len(result['relative_path']))
        col_widths[1] = max(col_widths[1], len(result['timestamp']))
        col_widths[2] = max(col_widths[2], len(result['label']))
        col_widths[3] = max(col_widths[3], len(result['input_fasta']))
        col_widths[4] = max(col_widths[4], len(result['results_file']))
        col_widths[5] = max(col_widths[5], len(result['total_reads']))
        col_widths[6] = max(col_widths[6], len(result['identical_matches']))
        col_widths[7] = max(col_widths[7], len(result['match_percent']))
    
    # Print table
    def print_separator():
        print('+' + '+'.join('-' * (width + 2) for width in col_widths) + '+')
    
    def print_row(values):
        row = '|'
        for i, value in enumerate(values):
            row += f' {value:<{col_widths[i]}} |'
        print(row)
    
    print_separator()
    print_row(headers)
    print_separator()
    
    for result in results:
        print_row([
            result['relative_path'],
            result['timestamp'],
            result['label'],
            result['input_fasta'],
            result['results_file'],
            result['total_reads'],
            result['identical_matches'],
            result['match_percent']
        ])
    
    print_separator()
    print(f"\nTotal files processed: {len(results)}")

def main():
    """Main function to process benchmark files and display results."""
    base_folder = "../data/reference"
    
    print(f"Processing benchmark files in: {base_folder}")
    print("=" * 50)
    
    results = process_benchmark_files(base_folder)
    
    if results:
        print("\nResults:")
        print_ascii_table(results)
    else:
        print("No valid benchmark files found.")

if __name__ == "__main__":
    main()

