#!/usr/bin/env python3

import os
import re
import glob
import csv
from pathlib import Path
from Bio import SeqIO
import pandas as pd
from collections import defaultdict
from multitax import SilvaTx

CSV_FILENAME = 'reference_evaluation.csv'
BASE_FOLDER = '../data/reference'

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

def load_multitax_silva():
    """Load MultiTax SILVA taxonomy system."""
    try:
        print("Loading SILVA accession-taxID mapping from taxmap file...")
        # Load taxmap file with primaryAccession and taxid columns
        taxmap_df = pd.read_csv(
            "taxmap_slv_ssu_ref_138.txt.gz", 
            sep='\t', 
            dtype={'taxid': str}
        )
        
        print(f"Loaded {len(taxmap_df)} total entries from taxmap file")
        
        # Remove duplicates by keeping first occurrence of each primaryAccession
        taxmap_df_unique = taxmap_df.drop_duplicates(subset=['primaryAccession'], keep='first')
        print(f"After removing duplicates: {len(taxmap_df_unique)} unique accessions")
        
        # Extract primaryAccession and taxid columns for MultiTax lookup
        acc_taxid = taxmap_df_unique[['primaryAccession', 'taxid']].set_index('primaryAccession')
        
        # Extract primaryAccession and organism_name for species-level information
        acc_organism = taxmap_df_unique[['primaryAccession', 'organism_name']].set_index('primaryAccession')
        
        print(f"Created {len(acc_taxid)} accession-taxID mappings")
        print(f"Created {len(acc_organism)} accession-organism mappings for species information")
        
        print("Initializing SILVA taxonomy system...")
        # Initialize with automatic SILVA download and lineage caching
        silva_tax = SilvaTx(
            files=["tax_slv_ssu_138.txt.gz"],
            build_node_children=True,
            build_name_nodes=True,
            build_rank_nodes=True,
            extended_names=True
        )
        
        print("SILVA taxonomy system initialized successfully")
        return acc_taxid, acc_organism, silva_tax
        
    except Exception as e:
        print(f"Error loading MultiTax SILVA: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None

def get_taxonomy_from_accession(accession, acc_taxid, acc_organism, silva_tax):
    """Get taxonomy lineage from accession using MultiTax."""
    try:
        # Strip everything past the first dot (exclude the dot as well)
        clean_accession = accession.split('.')[0] if '.' in accession else accession
        
        if clean_accession not in acc_taxid.index:
            return None
            
        taxid_result = acc_taxid.loc[clean_accession, 'taxid']
        
        # Handle case where we might get a Series (multiple matches despite deduplication)
        if isinstance(taxid_result, pd.Series):
            taxid = taxid_result.iloc[0]  # Take the first value
            print(f"Warning: Multiple taxids found for {clean_accession}, using first: {taxid}")
        else:
            taxid = taxid_result
        
        # Get lineage and ranks from MultiTax (up to genus level)
        lineage = silva_tax.name_lineage(taxid)
        ranks = silva_tax.rank_lineage(taxid)
        
        # Convert to dictionary mapping rank to name
        taxonomy_dict = {}
        for i, (name, rank) in enumerate(zip(lineage, ranks)):
            if rank and name:
                # Normalize rank names to match our standard levels
                rank_lower = rank.lower()
                if rank_lower in ['domain', 'phylum', 'class', 'order', 'family', 'genus']:
                    taxonomy_dict[rank_lower] = name
        
        # Add species information from organism_name if available
        if clean_accession in acc_organism.index:
            organism_result = acc_organism.loc[clean_accession, 'organism_name']
            
            # Handle case where we might get a Series
            if isinstance(organism_result, pd.Series):
                organism_name = organism_result.iloc[0]  # Take the first value
                print(f"Warning: Multiple organism names found for {clean_accession}, using first: {organism_name}")
            else:
                organism_name = organism_result
                
            if organism_name and pd.notna(organism_name):
                # Clean up organism name (remove parenthetical information)
                species_name = organism_name.split('(')[0].strip()
                if species_name:
                    taxonomy_dict['species'] = species_name
        
        return taxonomy_dict
        
    except Exception as e:
        print(f"Error getting taxonomy from accession {accession} (clean: {clean_accession}): {e}")
        return None

def analyze_single_record_identity(record, read_id, id_mapping):
    """Analyze identity for a single FASTA record. Returns (full_match, partial_match)."""
    if read_id not in id_mapping:
        return False, False
    
    mapping = id_mapping[read_id].lower()
    header_lower = record.description.lower()
    
    # Check for full mapping presence
    full_match = mapping in header_lower
    
    # Check for partial mapping (part before dot) presence
    partial_match = False
    if '.' in mapping:
        partial_mapping = mapping.split('.')[0]
        partial_match = partial_mapping in header_lower
    else:
        # If no dot, partial match is same as full match
        partial_match = full_match
    
    return full_match, partial_match

def analyze_single_record_taxonomy(record, read_id, id_mapping, acc_taxid, acc_organism, silva_tax):
    """Analyze taxonomic accuracy for a single FASTA record using MultiTax."""
    metrics_data = {}
    standard_levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    
    # Initialize metrics data for per-sequence tracking
    for level in standard_levels:
        metrics_data[level] = {
            'correct': 0,      # True Positives: correctly classified
            'incorrect': 0,    # False Positives: incorrectly classified  
            'not_classified': 0 # False Negatives: not classified
        }
    
    # Initialize fuzzy species matching metrics
    metrics_data['species_fuzzy'] = {
        'correct': 0,      # True Positives: correctly classified (fuzzy)
        'incorrect': 0,    # False Positives: incorrectly classified (fuzzy)
        'not_classified': 0 # False Negatives: not classified (fuzzy)
    }
    
    # Extract reference accession from header
    header = record.description
    reference_match = re.search(r'reference=([^ ,]+)', header)
    if not reference_match:
        # If we can't extract reference, count as not classified for all levels
        for level in standard_levels:
            metrics_data[level]['not_classified'] = 1
        metrics_data['species_fuzzy']['not_classified'] = 1
        return metrics_data
    
    reference_accession = reference_match.group(1)
    
    # Get true taxonomy from reference accession
    true_taxonomy = get_taxonomy_from_accession(reference_accession, acc_taxid, acc_organism, silva_tax)
    if not true_taxonomy:
        # If reference taxonomy not found, count as not classified for all levels
        for level in standard_levels:
            metrics_data[level]['not_classified'] = 1
        metrics_data['species_fuzzy']['not_classified'] = 1
        return metrics_data
    
    # Get predicted taxonomy (if any)
    predicted_taxonomy = None
    if read_id in id_mapping:
        predicted_accession = id_mapping[read_id]
        predicted_taxonomy = get_taxonomy_from_accession(predicted_accession, acc_taxid, acc_organism, silva_tax)
    
    # For each taxonomic level
    for level in standard_levels:
        # Get true taxon at this level
        true_taxon_at_level = true_taxonomy.get(level)
        
        if not true_taxon_at_level:
            # No true taxonomy at this level, skip
            continue
            
        # Get predicted taxon at this level (if we have a prediction)
        predicted_at_level = None
        if predicted_taxonomy:
            predicted_at_level = predicted_taxonomy.get(level)
        
        # Classify the prediction result
        if predicted_at_level:
            if predicted_at_level.lower().strip() == true_taxon_at_level.lower().strip():
                metrics_data[level]['correct'] = 1  # True Positive
            else:
                metrics_data[level]['incorrect'] = 1  # False Positive
        else:
            metrics_data[level]['not_classified'] = 1  # False Negative
    
    # Handle fuzzy species matching separately
    true_species = true_taxonomy.get('species')
    if true_species:
        predicted_species = None
        if predicted_taxonomy:
            predicted_species = predicted_taxonomy.get('species')
        
        if predicted_species:
            # Check for exact match first
            if predicted_species.lower().strip() == true_species.lower().strip():
                metrics_data['species_fuzzy']['correct'] = 1  # True Positive (exact)
            else:
                # Check for fuzzy match (one starts with the other)
                true_clean = true_species.lower().strip()
                pred_clean = predicted_species.lower().strip()
                
                # Only do fuzzy matching if both have length > 0
                if len(true_clean) > 0 and len(pred_clean) > 0:
                    if true_clean.startswith(pred_clean) or pred_clean.startswith(true_clean):
                        metrics_data['species_fuzzy']['correct'] = 1  # True Positive (fuzzy)
                    else:
                        metrics_data['species_fuzzy']['incorrect'] = 1  # False Positive
                else:
                    metrics_data['species_fuzzy']['incorrect'] = 1  # False Positive
        else:
            metrics_data['species_fuzzy']['not_classified'] = 1  # False Negative
    # If no true species, fuzzy matching is not applicable, leave at 0
    
    return metrics_data

def analyze_fasta(fasta_file, id_mapping, acc_taxid, acc_organism, silva_tax, label=""):
    """
    Analyze FASTA file once for both identity and taxonomic accuracy.
    Returns identity metrics and taxonomic accuracy metrics.
    """
    try:
        if not os.path.exists(fasta_file):
            print(f"Warning: FASTA file {fasta_file} not found")
            return 0, 0, 0, {}, {}
        
        # Initialize counters
        full_identity_matches = 0
        partial_identity_matches = 0
        total_reads = 0
        
        # Initialize taxonomic metrics aggregation for per-sequence tracking
        standard_levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        aggregated_taxonomy_data = {}
        for level in standard_levels:
            aggregated_taxonomy_data[level] = {
                'correct': 0,      # Total True Positives
                'incorrect': 0,    # Total False Positives
                'not_classified': 0 # Total False Negatives
            }
        
        # Initialize fuzzy species matching aggregation
        aggregated_taxonomy_data['species_fuzzy'] = {
            'correct': 0,      # Total True Positives (fuzzy)
            'incorrect': 0,    # Total False Positives (fuzzy)
            'not_classified': 0 # Total False Negatives (fuzzy)
        }
        
        # Process each record in the FASTA file
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description
            
            # Extract ID (number after > followed by space)
            id_match = re.match(r'^(\d+)\s', header)
            if not id_match:
                continue
            
            read_id = id_match.group(1)
            total_reads += 1
            
            # Analyze identity for this record
            full_match, partial_match = analyze_single_record_identity(record, read_id, id_mapping)
            if full_match:
                full_identity_matches += 1
            if partial_match:
                partial_identity_matches += 1
            
            # Analyze taxonomic accuracy for this record
            taxonomy_data = analyze_single_record_taxonomy(record, read_id, id_mapping, acc_taxid, acc_organism, silva_tax)
            
            # Aggregate taxonomic data
            for level in standard_levels:
                aggregated_taxonomy_data[level]['correct'] += taxonomy_data[level]['correct']
                aggregated_taxonomy_data[level]['incorrect'] += taxonomy_data[level]['incorrect']
                aggregated_taxonomy_data[level]['not_classified'] += taxonomy_data[level]['not_classified']
            
            # Aggregate fuzzy species data
            aggregated_taxonomy_data['species_fuzzy']['correct'] += taxonomy_data['species_fuzzy']['correct']
            aggregated_taxonomy_data['species_fuzzy']['incorrect'] += taxonomy_data['species_fuzzy']['incorrect']
            aggregated_taxonomy_data['species_fuzzy']['not_classified'] += taxonomy_data['species_fuzzy']['not_classified']
        
        print(f"Summary for {label} - {fasta_file}:")

        # Calculate final taxonomic metrics
        taxonomic_metrics = {}
        taxonomic_raw_counts = {}
        
        # Process standard levels
        for level in standard_levels:
            correct = aggregated_taxonomy_data[level]['correct']
            incorrect = aggregated_taxonomy_data[level]['incorrect']
            not_classified = aggregated_taxonomy_data[level]['not_classified']
            
            # Calculate precision, recall, F1 using per-sequence counts
            precision, recall, f1_score = calculate_sequence_precision_recall_f1(correct, incorrect, not_classified)
            
            taxonomic_metrics[level] = {
                'precision': precision,
                'recall': recall,
                'f1_score': f1_score
            }
            
            # Store raw counts
            taxonomic_raw_counts[level] = {
                'tp': correct,        # True Positives
                'fp': incorrect,      # False Positives  
                'fn': not_classified  # False Negatives
            }
            
            total_sequences = correct + incorrect + not_classified
            print(f"Level {level}: {correct} correct, {incorrect} incorrect, {not_classified} not classified (out of {total_sequences} total)")
        
        # Process fuzzy species matching
        correct = aggregated_taxonomy_data['species_fuzzy']['correct']
        incorrect = aggregated_taxonomy_data['species_fuzzy']['incorrect']
        not_classified = aggregated_taxonomy_data['species_fuzzy']['not_classified']
        
        precision, recall, f1_score = calculate_sequence_precision_recall_f1(correct, incorrect, not_classified)
        
        taxonomic_metrics['species_fuzzy'] = {
            'precision': precision,
            'recall': recall,
            'f1_score': f1_score
        }
        
        taxonomic_raw_counts['species_fuzzy'] = {
            'tp': correct,
            'fp': incorrect,
            'fn': not_classified
        }
        
        total_sequences = correct + incorrect + not_classified
        print(f"Level species_fuzzy: {correct} correct, {incorrect} incorrect, {not_classified} not classified (out of {total_sequences} total)")
        
        print(f"Identity analysis: {full_identity_matches} full matches, {partial_identity_matches} partial matches out of {total_reads} total reads")
        
        return full_identity_matches, partial_identity_matches, total_reads, taxonomic_metrics, taxonomic_raw_counts
    
    except Exception as e:
        print(f"Error analyzing FASTA file: {e}")
        import traceback
        traceback.print_exc()
        
        # Return default values
        taxonomic_metrics = {}
        taxonomic_raw_counts = {}
        for level in standard_levels:
            taxonomic_metrics[level] = {'precision': 0.0, 'recall': 0.0, 'f1_score': 0.0}
            taxonomic_raw_counts[level] = {'tp': 0, 'fp': 0, 'fn': 0}
        # Add species_fuzzy default values
        taxonomic_metrics['species_fuzzy'] = {'precision': 0.0, 'recall': 0.0, 'f1_score': 0.0}
        taxonomic_raw_counts['species_fuzzy'] = {'tp': 0, 'fp': 0, 'fn': 0}
        return 0, 0, 0, taxonomic_metrics, taxonomic_raw_counts

def calculate_sequence_precision_recall_f1(correct, incorrect, not_classified):
    """Calculate precision, recall, and F1 score using per-sequence counts."""
    if correct == 0 and incorrect == 0 and not_classified == 0:
        return 1.0, 1.0, 1.0  # Perfect if no sequences to classify
    
    if correct == 0:
        return 0.0, 0.0, 0.0  # No correct predictions
    
    # TP = correct, FP = incorrect, FN = not_classified
    precision = correct / (correct + incorrect) if (correct + incorrect) > 0 else 0.0
    recall = correct / (correct + not_classified) if (correct + not_classified) > 0 else 0.0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    
    return precision, recall, f1_score

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
            output_file = report_match.group(1)
            # Convert output.txt to headers.txt
            results_file = output_file.replace('output.txt', 'headers.txt')
    
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
        is_csv = label.startswith('sina') or label.startswith('kraken2')
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

def get_column_indices(label):
    """Get column indices for different command types."""
    if label.startswith('sina'):
        return 0, 3  # columns 1 and 4 (0-indexed)
    elif label.startswith('kraken2'):
        return 0, 2  # columns 2 and 5 (0-indexed)
    elif label.startswith('minimap2'):
        return 0, 5  # columns 1 and 6 (0-indexed)
    elif label.startswith('mothur'):
        return 0, 2  # columns 1 and 3 (0-indexed)
    else:
        return None, None

def process_benchmark_files(base_folder):
    """Process all benchmark files in the base folder recursively."""
    base_path = Path(base_folder)
    
    # Load existing processed entries
    processed_entries = load_existing_csv()
    
    # Load mothur to silva mappings
    seed_mapping, nr99_mapping = load_mothur_silva_mappings()
    
    # Load MultiTax SILVA taxonomy system
    acc_taxid, acc_organism, silva_tax = load_multitax_silva()
    
    # Check if MultiTax loaded successfully
    if acc_taxid is None or acc_organism is None or silva_tax is None:
        print("Warning: MultiTax SILVA system failed to load. Taxonomy analysis will be skipped.")
        print("Only identity analysis will be performed.")
    
    processed_count = 0
    skipped_count = 0
    
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
            
            # Get relative path from base folder (moved earlier)
            relative_path = str(log_file.parent.relative_to(base_path))
            
            # Check if this entry has already been processed
            entry_key = (relative_path, label, input_fasta or 'N/A')
            if entry_key in processed_entries:
                print(f"Skipping already processed entry: {entry_key}")
                skipped_count += 1
                continue
            
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
            
            if label == 'kraken2' and input_fasta == '10k':
                pass
            else:
                #continue
                pass

            if 'ncbi' in str(results_path):
                continue
            
            # Parse results file and get ID->mapping pairs
            col1_idx, col2_idx = get_column_indices(label)
            if col1_idx is None or col2_idx is None:
                print(f"Error: Unknown column indices for label {label}")
                continue
            
            print(f"Processing entry: {entry_key}")
            id_mapping = parse_results_file(results_path, label, col1_idx, col2_idx, seed_mapping, nr99_mapping)
            
            # Analyze FASTA file for identity and taxonomic accuracy
            total_reads = 0
            identical_matches = 0
            match_percent = "0.00"
            partial_match_percent = "0.00"
            
            # Calculate precision, recall, F1 metrics
            taxonomic_metrics = {}
            taxonomic_raw_counts = {}
            
            if input_fasta:
                fasta_file = f"{input_fasta}.fasta"
                
                # Only perform taxonomy analysis if MultiTax loaded successfully
                if acc_taxid is not None and acc_organism is not None and silva_tax is not None:
                    full_matches, partial_identity_matches, total, taxonomic_metrics, taxonomic_raw_counts = analyze_fasta(fasta_file, id_mapping, acc_taxid, acc_organism, silva_tax, label)
                else:
                    # Only perform identity analysis
                    full_matches, partial_identity_matches, total = analyze_fasta_identity_only(fasta_file, id_mapping)
                    # Initialize empty taxonomic metrics
                    for level in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                        taxonomic_metrics[level] = {'precision': 0.0, 'recall': 0.0, 'f1_score': 0.0}
                        taxonomic_raw_counts[level] = {'tp': 0, 'fp': 0, 'fn': 0}
                    # Initialize species_fuzzy for fallback case
                    taxonomic_metrics['species_fuzzy'] = {'precision': 0.0, 'recall': 0.0, 'f1_score': 0.0}
                    taxonomic_raw_counts['species_fuzzy'] = {'tp': 0, 'fp': 0, 'fn': 0}
                
                total_reads = total
                identical_matches = full_matches
                
                if total > 0:
                    percentage = (full_matches / total) * 100
                    match_percent = f"{percentage:.2f}"
                    
                    partial_percentage = (partial_identity_matches / total) * 100
                    partial_match_percent = f"{partial_percentage:.2f}"
            
            result_data = {
                'relative_path': relative_path,
                'timestamp': timestamp,
                'label': label,
                'input_fasta': input_fasta or 'N/A',
                'results_file': results_file,
                'total_reads': str(total_reads),
                'identical_matches': str(identical_matches),
                'match_percent': match_percent + '%',
                'partial_matches': str(partial_identity_matches),
                'partial_match_percent': partial_match_percent + '%'
            }
            
            # Add taxonomic metrics and raw counts to result data
            for level in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'species_fuzzy']:
                if level in taxonomic_metrics:
                    result_data[f'{level}_precision'] = f"{taxonomic_metrics[level]['precision']:.3f}"
                    result_data[f'{level}_recall'] = f"{taxonomic_metrics[level]['recall']:.3f}"
                    result_data[f'{level}_f1_score'] = f"{taxonomic_metrics[level]['f1_score']:.3f}"
                    result_data[f'{level}_tp'] = str(taxonomic_raw_counts[level]['tp'])
                    result_data[f'{level}_fp'] = str(taxonomic_raw_counts[level]['fp'])
                    result_data[f'{level}_fn'] = str(taxonomic_raw_counts[level]['fn'])
                else:
                    result_data[f'{level}_precision'] = "0.000"
                    result_data[f'{level}_recall'] = "0.000"
                    result_data[f'{level}_f1_score'] = "0.000"
                    result_data[f'{level}_tp'] = "0"
                    result_data[f'{level}_fp'] = "0"
                    result_data[f'{level}_fn'] = "0"
            
            # Append result to CSV immediately
            append_result_to_csv(result_data)
            
            # Add to processed entries to avoid reprocessing in current session
            processed_entries.add(entry_key)
            processed_count += 1
            
        except Exception as e:
            print(f"Error processing {log_file}: {e}")
    
    print(f"\nProcessing completed:")
    print(f"- Processed: {processed_count} new entries")
    print(f"- Skipped: {skipped_count} existing entries")


def analyze_fasta_identity_only(fasta_file, id_mapping):
    """
    Analyze FASTA file for identity only (when taxonomy system is not available).
    Returns only identity metrics.
    """
    try:
        if not os.path.exists(fasta_file):
            print(f"Warning: FASTA file {fasta_file} not found")
            return 0, 0, 0
        
        # Initialize counters
        full_identity_matches = 0
        partial_identity_matches = 0
        total_reads = 0
        
        # Process each record in the FASTA file
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description
            
            # Extract ID (number after > followed by space)
            id_match = re.match(r'^(\d+)\s', header)
            if not id_match:
                continue
            
            read_id = id_match.group(1)
            total_reads += 1
            
            # Analyze identity for this record
            full_match, partial_match = analyze_single_record_identity(record, read_id, id_mapping)
            if full_match:
                full_identity_matches += 1
            if partial_match:
                partial_identity_matches += 1
        
        print(f"Summary for {fasta_file} (identity only):")
        print(f"Identity analysis: {full_identity_matches} full matches, {partial_identity_matches} partial matches out of {total_reads} total reads")
        
        return full_identity_matches, partial_identity_matches, total_reads
    
    except Exception as e:
        print(f"Error analyzing FASTA file: {e}")
        import traceback
        traceback.print_exc()
        return 0, 0, 0

def save_results_to_csv(results, filename="reference_evaluation.csv"):
    """Save results to CSV file."""
    if not results:
        print("No results to save.")
        return
    
    try:
        # Create DataFrame from results
        df = pd.DataFrame(results)
        
        # Reorder columns to have main metrics first, then taxonomic metrics
        main_cols = ['relative_path', 'timestamp', 'label', 'input_fasta', 'results_file', 
                    'total_reads', 'identical_matches', 'match_percent', 'partial_matches', 'partial_match_percent']
        
        # Add taxonomic metric columns
        taxonomic_levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'species_fuzzy']
        taxonomic_cols = []
        for level in taxonomic_levels:
            taxonomic_cols.extend([f'{level}_precision', f'{level}_recall', f'{level}_f1_score', f'{level}_tp', f'{level}_fp', f'{level}_fn'])
        
        # Reorder columns
        all_cols = main_cols + taxonomic_cols
        existing_cols = [col for col in all_cols if col in df.columns]
        df = df[existing_cols]
        
        # Save to CSV
        df.to_csv(filename, index=False)
        print(f"Results saved to {filename}")
        
    except Exception as e:
        print(f"Error saving results to CSV: {e}")


def main():
    """Main function to process benchmark files and display results."""
    base_folder = BASE_FOLDER
    
    print(f"Processing benchmark files in: {base_folder}")
    print("=" * 50)
    
    process_benchmark_files(base_folder)
    
    print(f"\nResults are being saved incrementally to {CSV_FILENAME}")
    print("Processing complete!")

def load_existing_csv(filename=CSV_FILENAME):
    """Load existing CSV and return a set of processed entries (relative_path, label, input_fasta)."""
    processed_entries = set()
    
    if not os.path.exists(filename):
        return processed_entries
    
    try:
        df = pd.read_csv(filename)
        for _, row in df.iterrows():
            entry = (row['relative_path'], row['label'], row['input_fasta'])
            processed_entries.add(entry)
        print(f"Loaded {len(processed_entries)} existing entries from {filename}")
    except Exception as e:
        print(f"Error reading existing CSV {filename}: {e}")
    
    return processed_entries

def input_fasta_to_sort_key(input_fasta):
    """Convert input_fasta format (e.g., '10k', '5m') to sortable integer."""
    if not input_fasta or input_fasta == 'N/A':
        return 0  # Put N/A entries first
    
    # Extract number and suffix
    match = re.match(r'^(\d+)([km])$', input_fasta.lower())
    if not match:
        return 0  # Fallback for unexpected formats
    
    number = int(match.group(1))
    suffix = match.group(2)
    
    if suffix == 'k':
        return number * 1000
    elif suffix == 'm':
        return number * 1000000
    else:
        return number

def append_result_to_csv(result_data, filename=CSV_FILENAME):
    """Append a single result to the CSV file and maintain sorted order."""
    try:
        # Define the column order
        main_cols = ['relative_path', 'timestamp', 'label', 'input_fasta', 'results_file', 
                    'total_reads', 'identical_matches', 'match_percent', 'partial_matches', 'partial_match_percent']
        
        taxonomic_levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'species_fuzzy']
        taxonomic_cols = []
        for level in taxonomic_levels:
            taxonomic_cols.extend([f'{level}_precision', f'{level}_recall', f'{level}_f1_score', f'{level}_tp', f'{level}_fp', f'{level}_fn'])
        
        all_cols = main_cols + taxonomic_cols
        
        # Read existing data if file exists
        existing_data = []
        file_exists = os.path.exists(filename)
        
        if file_exists:
            try:
                df = pd.read_csv(filename)
                existing_data = df.to_dict('records')
            except Exception as e:
                print(f"Warning: Could not read existing CSV for sorting: {e}")
        
        # Add new result
        existing_data.append(result_data)
        
        # Sort by input_fasta (converted to numeric) then by label
        existing_data.sort(key=lambda x: (input_fasta_to_sort_key(x['input_fasta']), x['label']))
        
        # Write sorted data back to CSV
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=all_cols)
            writer.writeheader()
            writer.writerows(existing_data)
            
        print(f"Appended and sorted result in {filename}")
        
    except Exception as e:
        print(f"Error appending result to CSV: {e}")

if __name__ == "__main__":
    main()

