#!/usr/bin/env python3

import sys
import multiprocessing
import argparse
import time
import os
import csv
from Bio import Align, SeqIO, AlignIO
from Bio.Align import AlignInfo
from Bio.motifs import Motif
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
from functools import partial


COLUMN_CHUNK_SIZE = 100
ROW_CHUNK_SIZE = 100

# Debug settings
DEBUG_TIMING = True  # Set to False to disable timing measurements
DEBUG_PROGRESS = False  # Set to False to disable progress reporting

# RNA scoring matrix for alignment evaluation. Adapted from the NUC44 scoring matrix.
# The gap penalty is 0, as it is not used in the scoring. The gaps are scored separately, based on opens and extensions.
RNA_SCORING_MATRIX = {
    'A': {'A': 5, 'C': -4, 'G': -4, 'U': -4, '-': 0},
    'C': {'A': -4, 'C': 5, 'G': -4, 'U': -4, '-': 0},
    'G': {'A': -4, 'C': -4, 'G': 5, 'U': -4, '-': 0},
    'U': {'A': -4, 'C': -4, 'G': -4, 'U': 5, '-': 0},
    '-': {'A': 0, 'C': 0, 'G': 0, 'U': 0, '-': 0}
}

# RNA scoring matrix for alignment evaluation with built-in gap penalty.
RNA_SCORING_MATRIX_WITH_GAP_PENALTY = {
    'A': {'A': 5, 'C': -4, 'G': -4, 'U': -4, '-': -10},
    'C': {'A': -4, 'C': 5, 'G': -4, 'U': -4, '-': -10},
    'G': {'A': -4, 'C': -4, 'G': 5, 'U': -4, '-': -10},
    'U': {'A': -4, 'C': -4, 'G': -4, 'U': 5, '-': -10},
    '-': {'A': -10, 'C': -10, 'G': -10, 'U': -10, '-': 0}
}

GAP_OPEN_PENALTY = -10
GAP_EXTENSION_PENALTY = -0.5

# Weights for the composite score
ALPHA = 1.0  # Weight for standard score per pair
BETA = 1.0   # Weight for gap penalty per gap

def read_alignment(alignment_file):
    """
    Reads a multiple sequence alignment file.
    
    Args:
        alignment_file (str): Path to the alignment file in FASTA format
        
    Returns:
        Bio.Align.MultipleSeqAlignment: The alignment object
    """
    try:
        alignment = Align.read(alignment_file, "fasta")
        return alignment
    except Exception as e:
        print(f"Error reading alignment file: {e}", file=sys.stderr)
        sys.exit(1)

def process_columns_chunk(chunk_info):
    """
    Process a chunk of columns to count nucleotides and gaps.
    
    Args:
        chunk_info (tuple): Tuple containing (start_col, end_col, alignment_file)
        
    Returns:
        list: List of tuples (col_idx, Counter)
    """
    start_col, end_col, alignment_file = chunk_info
    results = []
    
    try:
        # Read the sequences from file
        sequences = []
        with open(alignment_file, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            sequences = [str(record.seq) for record in records]
        
        # Process each column in this chunk
        for col_idx in range(start_col, end_col):
            if DEBUG_PROGRESS and col_idx % 100 == 0:
                print(f"Processing column {col_idx}.")
                
            # Extract the column by accessing the corresponding position in each sequence
            column = [seq[col_idx] if col_idx < len(seq) else '-' for seq in sequences]
            counts = Counter(column)
            
            # Get counts for each nucleotide (case-insensitive)
            a_count = counts.get('A', 0) + counts.get('a', 0)
            c_count = counts.get('C', 0) + counts.get('c', 0)
            g_count = counts.get('G', 0) + counts.get('g', 0)
            # Also process T and t as U, just in case the input is not 100% RNA
            u_count = counts.get('U', 0) + counts.get('u', 0) + counts.get('T', 0) + counts.get('t', 0)
            gap_count = counts.get('-', 0)
            
            # Create a new Counter with normalized counts
            normalized_counts = Counter({'A': a_count, 'C': c_count, 'G': g_count, 'U': u_count, '-': gap_count})
            results.append((col_idx, normalized_counts))
            
    except Exception as e:
        print(f"Error processing columns {start_col}-{end_col}: {str(e)}")
        
    return results

def get_column_counts(alignment, alignment_file, num_threads):
    """
    Counts nucleotides and gaps in each column of the alignment in parallel.
    
    Args:
        alignment: Bio.Align.MultipleSeqAlignment object (only used to get dimensions)
        alignment_file: Path to the original alignment file
        num_threads (int): Number of threads to use for parallel processing
        
    Returns:
        list: List of Counter objects with counts for each column
    """
    total_columns = alignment.length
    
    # Create chunks of columns
    chunk_size = COLUMN_CHUNK_SIZE  # Process X columns at a time
    chunks = []
    
    for i in range(0, total_columns, chunk_size):
        end_idx = min(i + chunk_size, total_columns)
        chunks.append((i, end_idx, alignment_file))
    
    # Use configured number of threads
    effective_threads = min(num_threads, os.cpu_count())
    
    with multiprocessing.Pool(processes=effective_threads) as pool:
        chunk_results = pool.map(process_columns_chunk, chunks)
    
    # Flatten and sort the results
    flattened = [item for sublist in chunk_results for item in sublist]
    flattened.sort(key=lambda x: x[0])  # Sort by column index
    
    # Extract just the counters in order
    column_counts = [counts for _, counts in flattened]
    
    return column_counts

def calculate_column_score(counts, scoring_matrix=RNA_SCORING_MATRIX):
    """
    Calculate sum-of-pairs score for a column.
    
    Args:
        counts (Counter): Counts of nucleotides in the column
        scoring_matrix (dict): Scoring matrix to use for calculation
        
    Returns:
        int: The score for the column
    """
    column_score = 0
    nucleotides = {'A': counts['A'], 'C': counts['C'], 'G': counts['G'], 'U': counts['U'], '-': counts['-']}
    
    # Compare each pair of nucleotides
    for n1 in nucleotides:
        for n2 in nucleotides:
            if n1 <= n2:  # Avoid counting pairs twice
                if n1 == n2:
                    # For same nucleotide pairs, use (count * (count-1))/2 to get number of pairs
                    count = nucleotides[n1]
                    pair_count = (count * (count - 1)) // 2
                    column_score += scoring_matrix[n1][n2] * pair_count
                else:
                    # For different nucleotides, multiply the counts
                    column_score += scoring_matrix[n1][n2] * nucleotides[n1] * nucleotides[n2]
    
    return column_score

def calculate_alignment_scores(column_counts):
    """
    Calculate alignment scores using both scoring matrices.
    
    Args:
        column_counts (list): List of Counter objects with counts for each column
        
    Returns:
        tuple: (standard_score, gap_penalty_score, total_nucleotides, total_pairs, total_gaps) 
    """
    standard_sum_score = 0
    gap_penalty_sum_score = 0
    total_nucleotides = 0
    total_pairs = 0
    total_gaps = 0
    
    for col_idx, counts in enumerate(column_counts):
        # Calculate both scores
        standard_score = calculate_column_score(counts, RNA_SCORING_MATRIX)
        gap_penalty_score = calculate_column_score(counts, RNA_SCORING_MATRIX_WITH_GAP_PENALTY)
        
        # Add to totals
        standard_sum_score += standard_score
        gap_penalty_sum_score += gap_penalty_score
        
        # Count nucleotides and gaps
        nucleotides_in_column = counts['A'] + counts['C'] + counts['G'] + counts['U']
        gaps_in_column = counts['-']
        total_nucleotides += nucleotides_in_column
        total_gaps += gaps_in_column
        
        # Count pairs in this column
        total_chars = nucleotides_in_column + gaps_in_column
        # Number of pairs = n(n-1)/2 where n is total characters
        pairs_in_column = (total_chars * (total_chars - 1)) // 2
        total_pairs += pairs_in_column
    
    return (standard_sum_score, gap_penalty_sum_score, total_nucleotides, total_pairs, total_gaps)

def process_rows_chunk(chunk_info):
    """
    Process a chunk of rows to count gap opens and extensions.
    
    Args:
        chunk_info (tuple): Tuple containing (start_idx, end_idx, alignment_file)
        
    Returns:
        list: List of tuples (row_index, gap_opens, gap_extensions)
    """
    start_idx, end_idx, alignment_file = chunk_info
    results = []
    
    try:
        # Read only the needed sequences from the file
        # This avoids passing large alignment objects between processes
        alignment = None
        sequences = []
        
        # Load just this subset of sequences
        with open(alignment_file, "r") as handle:
            all_records = list(SeqIO.parse(handle, "fasta"))
            sequences = all_records[start_idx:end_idx]
        
        for i, record in enumerate(sequences):
            row_idx = start_idx + i
            if DEBUG_PROGRESS and row_idx % 100 == 0:
                print(f"Processing row {row_idx}")
            
            # Get sequence string from the record
            sequence = str(record.seq)
            
            gap_opens = 0
            gap_extensions = 0
            in_gap = False
            
            # Iterate through the sequence
            for char in sequence:
                if char == '-':
                    if not in_gap:
                        # This is a new gap opening
                        gap_opens += 1
                        in_gap = True
                    else:
                        # This is a gap extension
                        gap_extensions += 1
                else:
                    # Not a gap
                    in_gap = False
            
            results.append((row_idx, gap_opens, gap_extensions))
            
    except Exception as e:
        print(f"Error processing chunk {start_idx}-{end_idx}: {str(e)}")
    
    return results

def get_row_gap_stats(alignment, alignment_file, num_threads):
    """
    Analyze rows of the alignment to count gap opens and extensions in parallel.
    
    Args:
        alignment: Bio.Align.MultipleSeqAlignment object (only used to get dimensions)
        alignment_file: Path to the original alignment file
        num_threads (int): Number of threads to use for parallel processing
        
    Returns:
        list: List of tuples (row_index, gap_opens, gap_extensions)
    """
    total_rows = len(alignment)
    
    # Create chunks of rows, more efficiently using file-based loading
    chunk_size = ROW_CHUNK_SIZE  # Process X sequences at a time
    chunks = []
    
    for i in range(0, total_rows, chunk_size):
        end_idx = min(i + chunk_size, total_rows)
        chunks.append((i, end_idx, alignment_file))
    
    # Use configured number of threads
    effective_threads = min(num_threads, os.cpu_count())
    
    with multiprocessing.Pool(processes=effective_threads) as pool:
        chunk_results = pool.map(process_rows_chunk, chunks)
    
    # Flatten the results
    results = [item for sublist in chunk_results for item in sublist]
    
    # Sort by row index for consistency
    results.sort(key=lambda x: x[0])
    
    return results

def calculate_consensus(column_counts, allow_gaps=True, threshold=0.7, ambiguous_character='n'):
    """
    Calculate consensus sequence with configurable gap handling.
    
    Args:
        column_counts: List of Counter objects with counts for each column
        allow_gaps: Whether to consider gaps in frequency calculations
        threshold: Minimum frequency threshold for inclusion (0-1)
        ambiguous_character: Character to use when no nucleotide meets threshold
        
    Returns:
        str: Consensus sequence
    """
    consensus = []
    
    for counts in column_counts:
        # Calculate total including gaps
        total = counts['A'] + counts['C'] + counts['G'] + counts['U'] + counts['-']

        # If there are no valid nucleotides in this column
        if total == 0:
            consensus.append('-')
            continue
            
        # Check each nucleotide against threshold
        consensus_base = ambiguous_character  # Default if no base meets threshold
        
        for base in ['A', 'C', 'G', 'U']:
            if counts[base] / total >= threshold:
                consensus_base = base
                break
                
        # Check gaps only if we're allowing them
        if allow_gaps and counts['-'] / total >= threshold:
            consensus_base = '-'
        
        consensus.append(consensus_base)
    
    return ''.join(consensus)

def save_column_stats(column_counts, output_file):
    """Save column statistics to a CSV file."""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Write header
        writer.writerow(['Column', 'A', 'C', 'G', 'U', 'Gaps'])
        
        # Write data for each column
        for i, counts in enumerate(column_counts):
            writer.writerow([
                i+1,
                counts['A'],
                counts['C'],
                counts['G'],
                counts['U'],
                counts['-']
            ])

def save_row_gap_stats(row_gap_stats, output_file):
    """Save row gap statistics to a CSV file."""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Write header
        writer.writerow(['Row', 'Gap Opens', 'Gap Extensions'])
        
        # Write data for each row
        for row_idx, gap_opens, gap_extensions in row_gap_stats:
            writer.writerow([row_idx+1, gap_opens, gap_extensions])

def save_consensus_to_fasta(consensus_sequence, output_file, description="Consensus Sequence"):
    """Save a consensus sequence to a FASTA file."""
    record = SeqRecord(Seq(consensus_sequence), id="consensus", description=description)
    with open(output_file, 'w') as f:
        f.write(f">{record.description}\n")
        
        # Write sequence in chunks of 80 characters
        for i in range(0, len(consensus_sequence), 80):
            f.write(consensus_sequence[i:i+80] + "\n")

def save_stats_to_csv(stats, output_file):
    """Save statistics to a CSV file."""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write each stat as a separate row with name and value
        for name, value in stats.items():
            writer.writerow([name, value])

def parse_args():
    """Parse command line arguments."""
    # Determine default number of threads (cores - 1, minimum 1)
    default_threads = max(1, multiprocessing.cpu_count() - 1)
    
    parser = argparse.ArgumentParser(description="Evaluate RNA sequence alignments")
    parser.add_argument('--input', '-i', required=True, help='Path to the alignment file in FASTA format')
    parser.add_argument('--threads', '-t', type=int, default=default_threads,
                        help=f"Number of threads to use (default: {default_threads})")
    parser.add_argument('--output-base-name', '--base-output-name', '-o', '--output', required=True, type=str, 
                        help="Base name for output files. If provided, stats will be saved to files")
    return parser.parse_args()

def calculate_tc_score(column_counts, threshold=0.7, allow_gaps=False):
    """
    Calculate the Total Column (TC) score, which is the fraction of columns 
    that have a clear majority nucleotide above the specified threshold.
    
    Args:
        column_counts: List of Counter objects with counts for each column
        threshold: Minimum frequency threshold for a nucleotide to be considered dominant (default: 0.7)
        
    Returns:
        float: TC score between 0 and 1
    """
    total_columns = len(column_counts)
    identical_columns = 0

    if allow_gaps:
        base_string = 'ACGU-'
    else:
        base_string = 'ACGU'

    # Check each column for a dominant nucleotide
    for counts in column_counts:
        total = sum(counts.values())
        
        # Find the majority nucleotide
        max_count = 0

        for base in base_string:
            if counts[base] > max_count:
                max_count = counts[base]
        
        # If the dominant base exceeds the threshold, count this column
        if max_count / total >= threshold:
            identical_columns += 1
    
    # Calculate the TC score
    tc_score = identical_columns / total_columns if total_columns > 0 else 0
    return tc_score

def calculate_column_consistency_score(column_counts):
    """
    Calculate the Column Consistency Score, which is the average of the ratio between
    the predominant base count and total bases for each column.
    
    Args:
        column_counts: List of Counter objects with counts for each column
        
    Returns:
        float: Average consistency score between 0 and 1
    """
    total_columns = len(column_counts)
    if total_columns == 0:
        return 0
        
    column_scores = []
    base_string = 'ACGU'  # Excluding gaps
    
    # Calculate score for each column
    for counts in column_counts:
        total = sum(counts.values())
        if total == 0:
            continue
            
        # Find the majority nucleotide
        max_count = 0
        for base in base_string:
            if counts[base] > max_count:
                max_count = counts[base]
        
        # Calculate this column's score
        column_score = max_count / total
        column_scores.append(column_score)
    
    # Calculate the average consistency score
    consistency_score = sum(column_scores) / len(column_scores) if column_scores else 0
    return consistency_score

def evaluate_alignment(alignment_file, num_threads, output_base_name=None):
    """
    Evaluates a multiple sequence alignment.
    
    Args:
        alignment_file (str): Path to the alignment file in FASTA format
        num_threads (int): Number of threads to use for parallel processing
        output_base_name (str, optional): Base name for output files
    """
    # Read the alignment (just to get dimensions)
    alignment = read_alignment(alignment_file)
    print(f"Evaluating alignment using {num_threads} threads...")
    print(f"Alignment has {len(alignment)} sequences and {alignment.length} columns")    
    
    # Get column counts in parallel and time it if debug timing is enabled
    column_time_start = time.time() if DEBUG_TIMING else 0
    column_counts = get_column_counts(alignment, alignment_file, num_threads)
    if DEBUG_TIMING:
        column_time_end = time.time()
        print(f"Column processing took {column_time_end - column_time_start:.2f} seconds")
    
    # Calculate both types of consensus using the unified function with 70% threshold
    with_gap_70 = calculate_consensus(column_counts, allow_gaps=True, threshold=0.7)
    no_gap_70 = calculate_consensus(column_counts, allow_gaps=False, threshold=0.7)
    
    # Calculate both types of consensus using the unified function with 50% threshold
    with_gap_50 = calculate_consensus(column_counts, allow_gaps=True, threshold=0.5)
    no_gap_50 = calculate_consensus(column_counts, allow_gaps=False, threshold=0.5)
    
    # Print consensus sequences
    print("\nNo-Gap Consensus (70% threshold, gaps ignored):")
    # Print in blocks of 80 characters for readability
    for i in range(0, len(no_gap_70), 80):
        print(no_gap_70[i:i+80])
    
    print("\nGap Consensus (70% threshold, gaps included):")
    for i in range(0, len(with_gap_70), 80):
        print(with_gap_70[i:i+80])
    
    print("\nNo-Gap Consensus (50% threshold, gaps ignored):")
    for i in range(0, len(no_gap_50), 80):
        print(no_gap_50[i:i+80])
    
    print("\nGap Consensus (50% threshold, gaps included):")
    for i in range(0, len(with_gap_50), 80):
        print(with_gap_50[i:i+80])

    print("\n")
    
    standard_score, gap_penalty_score, total_nucleotides, total_pairs, total_gaps = calculate_alignment_scores(column_counts)
    
    # Get row gap statistics in parallel and time it if debug timing is enabled
    row_time_start = time.time() if DEBUG_TIMING else 0
    row_gap_stats = get_row_gap_stats(alignment, alignment_file, num_threads)
    if DEBUG_TIMING:
        row_time_end = time.time()
        print(f"Row gap analysis took {row_time_end - row_time_start:.2f} seconds")
    
    # Calculate and display gap statistics
    total_gap_opens = sum(stats[1] for stats in row_gap_stats)
    total_gap_extensions = sum(stats[2] for stats in row_gap_stats)
    avg_gap_opens = total_gap_opens / len(alignment)
    avg_gap_extensions = total_gap_extensions / len(alignment)
    avg_gap_opens_per_nucleotide = total_gap_opens / total_nucleotides if total_nucleotides > 0 else 0
    avg_gap_extensions_per_nucleotide = total_gap_extensions / total_nucleotides if total_nucleotides > 0 else 0
    
    print("\nGap Statistics:")
    print("-" * 50)
    print(f"Total gap opens: {total_gap_opens}")
    print(f"Total gap extensions: {total_gap_extensions}")
    print(f"Average gap opens per sequence: {avg_gap_opens:.2f}")
    print(f"Average gap extensions per sequence: {avg_gap_extensions:.2f}")
    print(f"Average gap opens per nucleotide: {avg_gap_opens_per_nucleotide:.4f}")
    print(f"Average gap extensions per nucleotide: {avg_gap_extensions_per_nucleotide:.4f}")

    # Calculate TC score with different thresholds
    tc_score_10 = calculate_tc_score(column_counts, threshold=0.1)
    tc_score_20 = calculate_tc_score(column_counts, threshold=0.2)
    tc_score_30 = calculate_tc_score(column_counts, threshold=0.3)
    tc_score_40 = calculate_tc_score(column_counts, threshold=0.4)
    tc_score_50 = calculate_tc_score(column_counts, threshold=0.5)
    tc_score_60 = calculate_tc_score(column_counts, threshold=0.6)
    tc_score_70 = calculate_tc_score(column_counts, threshold=0.7)
    tc_score_80 = calculate_tc_score(column_counts, threshold=0.8)
    tc_score_90 = calculate_tc_score(column_counts, threshold=0.9)
    tc_score_100 = calculate_tc_score(column_counts, threshold=1.0)
    
    # Calculate column consistency score
    column_consistency_score = calculate_column_consistency_score(column_counts)
    
    # Summary statistics
    print("\nSummary Statistics:")
    print("-" * 50)
    print(f"Alignment dimensions: {len(alignment)} sequences × {alignment.length} columns")
    print(f"Total nucleotides (excluding gaps): {total_nucleotides}")
    print(f"Total gaps: {total_gaps}")
    print(f"Total character pairs: {total_pairs}")
    print(f"Gap Fraction: {total_gaps / (total_gaps + total_nucleotides):.2%}")
    print(f"TC score (threshold=10%): {tc_score_10:.4f}")
    print(f"TC score (threshold=20%): {tc_score_20:.4f}")
    print(f"TC score (threshold=30%): {tc_score_30:.4f}")
    print(f"TC score (threshold=40%): {tc_score_40:.4f}")
    print(f"TC score (threshold=50%): {tc_score_50:.4f}")
    print(f"TC score (threshold=60%): {tc_score_60:.4f}")
    print(f"TC score (threshold=70%): {tc_score_70:.4f}")
    print(f"TC score (threshold=80%): {tc_score_80:.4f}")
    print(f"TC score (threshold=90%): {tc_score_90:.4f}")
    print(f"TC score (threshold=100%): {tc_score_100:.4f}")
    print(f"Column consistency score: {column_consistency_score:.4f}")

    # Calculate normalization factors
    num_columns = alignment.length
    num_rows = len(alignment)
    
    total_gap_penalty = total_gap_opens * GAP_OPEN_PENALTY + total_gap_extensions * GAP_EXTENSION_PENALTY
    
    # Calculate gap penalty normalizations
    gap_penalty_per_sequence = total_gap_penalty / num_rows
    gap_penalty_per_nucleotide = total_gap_penalty / total_nucleotides if total_nucleotides > 0 else 0
    gap_penalty_per_pair = total_gap_penalty / total_pairs if total_pairs > 0 else 0
    gap_penalty_per_gap = total_gap_penalty / total_gaps if total_gaps > 0 else 0
    
    # Calculate standard score normalizations
    standard_score_per_column = standard_score / num_columns
    standard_score_per_nucleotide = standard_score / total_nucleotides if total_nucleotides > 0 else 0
    standard_score_per_pair = standard_score / total_pairs if total_pairs > 0 else 0
    
    # Calculate new combined scores
    combined_score_per_pair = standard_score_per_pair + gap_penalty_per_pair
    combined_score_per_nucleotide = standard_score_per_nucleotide + gap_penalty_per_nucleotide
    composite_score = (ALPHA * standard_score_per_pair) + (BETA * gap_penalty_per_gap)
    
    # Calculate and display combined scores
    print("\nStandard Score (matrix with no gap penalty):")
    print("-" * 50)
    print(f"Total: {standard_score}")
    print(f"  - Per column: {standard_score_per_column:.4f}")
    print(f"  - Per nucleotide: {standard_score_per_nucleotide:.4f}")
    print(f"  - Per pair: {standard_score_per_pair:.4f}")
    
    print("\nGap Penalty (based on opens/extensions):")
    print("-" * 50)
    print(f"Total: {total_gap_penalty}")
    print(f"  - Per sequence: {gap_penalty_per_sequence:.4f}")
    print(f"  - Per nucleotide: {gap_penalty_per_nucleotide:.4f}")
    print(f"  - Per pair: {gap_penalty_per_pair:.4f}")
    print(f"  - Per gap: {gap_penalty_per_gap:.4f}")
    
    print("\nCombined Scores:")
    print("-" * 50)
    print(f"Per pair (standard per pair + gap penalty per pair): {combined_score_per_pair:.4f}")
    print(f"Per nucleotide (standard per nucleotide + gap penalty per nucleotide): {combined_score_per_nucleotide:.4f}")
    print(f"Composite score (α={ALPHA} × standard per pair + β={BETA} × gap penalty per gap): {composite_score:.4f}")
    
    print("\nAlternative Score (matrix with built-in gap penalty):")
    print("-" * 50)
    print(f"Total: {gap_penalty_score}")
    print(f"  - Per column: {gap_penalty_score / num_columns:.4f}")
    print(f"  - Per nucleotide: {gap_penalty_score / total_nucleotides:.4f}")
    print(f"  - Per pair: {gap_penalty_score / total_pairs:.4f}")
    
    # Save files if output_base_name is provided
    if output_base_name:
        print(f"\nSaving output files with base name: {output_base_name}")
        
        # Save column statistics
        column_stats_file = f"{output_base_name}_column_stats.csv"
        save_column_stats(column_counts, column_stats_file)
        print(f"Saved column statistics to {column_stats_file}")
        
        # Save row gap statistics
        row_gap_stats_file = f"{output_base_name}_row_gap_stats.csv"
        save_row_gap_stats(row_gap_stats, row_gap_stats_file)
        print(f"Saved row gap statistics to {row_gap_stats_file}")
        
        # Save consensus sequences to FASTA files
        save_consensus_to_fasta(no_gap_70, f"{output_base_name}_consensus_70_no_gaps.fasta", 
                               "70% Threshold Consensus (gaps ignored)")
        save_consensus_to_fasta(with_gap_70, f"{output_base_name}_consensus_70_with_gaps.fasta", 
                               "70% Threshold Consensus (gaps included)")
        save_consensus_to_fasta(no_gap_50, f"{output_base_name}_consensus_50_no_gaps.fasta", 
                               "50% Threshold Consensus (gaps ignored)")
        save_consensus_to_fasta(with_gap_50, f"{output_base_name}_consensus_50_with_gaps.fasta", 
                               "50% Threshold Consensus (gaps included)")
        print(f"Saved consensus sequences to {output_base_name}_consensus_*.fasta files")
        
        # Save all statistics to a CSV file
        stats = {
            "Total Rows": num_rows,
            "Total Columns": num_columns,
            "Total Nucleotides": total_nucleotides,
            "Total Gaps": total_gaps,
            "Total Character Pairs": total_pairs,
            "Gap Fraction": total_gaps / (total_gaps + total_nucleotides),
            "Total Gap Opens": total_gap_opens,
            "Total Gap Extensions": total_gap_extensions,
            "Average Gap Opens Per Sequence": avg_gap_opens,
            "Average Gap Extensions Per Sequence": avg_gap_extensions,
            "Average Gap Opens Per Nucleotide": avg_gap_opens_per_nucleotide,
            "Average Gap Extensions Per Nucleotide": avg_gap_extensions_per_nucleotide,
            "Standard Score": standard_score,
            "Standard Score Per Column": standard_score_per_column,
            "Standard Score Per Nucleotide": standard_score_per_nucleotide,
            "Standard Score Per Pair": standard_score_per_pair,
            "Gap Penalty": total_gap_penalty,
            "Gap Penalty Per Sequence": gap_penalty_per_sequence,
            "Gap Penalty Per Nucleotide": gap_penalty_per_nucleotide,
            "Gap Penalty Per Pair": gap_penalty_per_pair,
            "Gap Penalty Per Gap": gap_penalty_per_gap,
            "Combined Score Per Pair (standard per pair + gap penalty per pair)": combined_score_per_pair,
            "Combined Score Per Nucleotide (standard per nucleotide + gap penalty per nucleotide)": combined_score_per_nucleotide,
            "Alternative Score (With Gap Penalty Matrix)": gap_penalty_score,
            "Alternative Score Per Column": gap_penalty_score / num_columns,
            "Alternative Score Per Nucleotide": gap_penalty_score / total_nucleotides,
            "Alternative Score Per Pair": gap_penalty_score / total_pairs,
            "TC Score (10% threshold)": tc_score_10,
            "TC Score (20% threshold)": tc_score_20,
            "TC Score (30% threshold)": tc_score_30,
            "TC Score (40% threshold)": tc_score_40,
            "TC Score (50% threshold)": tc_score_50,
            "TC Score (60% threshold)": tc_score_60,
            "TC Score (70% threshold)": tc_score_70,
            "TC Score (80% threshold)": tc_score_80,
            "TC Score (90% threshold)": tc_score_90,
            "TC Score (100% threshold)": tc_score_100,
            "Column Consistency Score": column_consistency_score
        }
        stats_file = f"{output_base_name}_stats.csv"
        save_stats_to_csv(stats, stats_file)
        print(f"Saved statistics to {stats_file}")

if __name__ == "__main__":
    args = parse_args()
    evaluate_alignment(args.input, args.threads, args.output_base_name)
