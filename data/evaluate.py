#! /usr/bin/env python3

import argparse
import time
import csv
import re
import os
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple, Optional
import numpy as np

from Bio import SeqIO
from Bio.Align import PairwiseAligner, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator
import Levenshtein
import fastdtw
import jellyfish
import pylcs

DEBUG_TIMING = True

def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate sequence alignments against reference and source data.')
    
    parser.add_argument('--input', '-i', required=True, help='Input file in FASTA format')
    parser.add_argument('--reference-alignment', '--reference', '-r', required=True, 
                        help='Reference alignment file in FASTA format')
    parser.add_argument('--source-data', '--source', '-s', required=True, 
                        help='Source data file in FASTA format')
    parser.add_argument('--output-base-name', '--base-output-name', '-o', '--output', required=True, 
                        help='Base name for output files')
    parser.add_argument('--threads', '-t', type=int, default=max(1, multiprocessing.cpu_count() - 1),
                        help='Number of threads to use for parallelization (default: number of CPU cores - 1)')
    parser.add_argument('--with-dtw', '--compute-dtw', '--dtw', action='store_true', default=False,
                        help='Enable DTW (Dynamic Time Warping) calculations (default: False)')
    
    return parser.parse_args()

def extract_reference_id(header: str) -> str:
    """Extract the first reference ID from a header string."""
    match = re.search(r'reference=([^,\s]+)', header)
    if match:
        return match.group(1)
    return None

def seq_to_numeric(seq: str) -> List[int]:
    """Convert a sequence to numeric representation for DTW."""
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 3, 'N': 4, '-': 5}
    return [mapping.get(c.upper(), 6) for c in seq]

def compute_metrics(seq1: str, seq2: str, with_dtw: bool = False, timers: Dict[str, float] = None) -> Dict[str, float]:
    """Compute all comparison metrics between two sequences."""
    metrics = {}
    
    # Convert sequences to uppercase
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    # Replace U with T to ensure DNA-compatible sequences
    seq1_dna = seq1.replace('U', 'T').replace('u', 't')
    seq2_dna = seq2.replace('U', 'T').replace('u', 't')
    
    # Levenshtein distance
    if timers is not None:
        start = time.time()
    try:
        lev_dist = Levenshtein.distance(seq1, seq2)
    except Exception as e:
        print(f"Error in Levenshtein.distance with seq1: {seq1}, seq2: {seq2}")
        raise
    metrics['levenshtein'] = lev_dist
    metrics['normalized_levenshtein'] = lev_dist / max(len(seq1), len(seq2))
    # Normalize against ref sequence length
    metrics['normalized_levenshtein_ref'] = lev_dist / len(seq2)
    if timers is not None:
        timers['levenshtein'] += time.time() - start
    
    # DTW (Dynamic Time Warping) - Only compute if requested
    if with_dtw:
        if timers is not None:
            start = time.time()
        try:
            seq1_num = seq_to_numeric(seq1)
            seq2_num = seq_to_numeric(seq2)
            dtw_dist, _ = fastdtw.fastdtw(seq1_num, seq2_num)
        except Exception as e:
            print(f"Error in fastdtw with seq1: {seq1}, seq2: {seq2}")
            raise
        metrics['dtw'] = dtw_dist
        metrics['normalized_dtw'] = dtw_dist / (len(seq1) + len(seq2))
        # Normalize against ref sequence length
        metrics['normalized_dtw_ref'] = dtw_dist / len(seq2)
        if timers is not None:
            timers['dtw'] += time.time() - start
    
    # Jaro-Winkler Similarity
    if timers is not None:
        start = time.time()
    try:
        jaro = jellyfish.jaro_winkler_similarity(seq1, seq2)
    except Exception as e:
        print(f"Error in jaro_winkler_similarity with seq1: {seq1}, seq2: {seq2}")
        raise
    metrics['jaro_winkler'] = jaro
    # Already normalized between 0 and 1
    if timers is not None:
        timers['jaro_winkler'] += time.time() - start
    
    # LCS (Longest Common Subsequence) sequence length
    if timers is not None:
        start = time.time()
    try:
        lcs = pylcs.lcs_string_length(seq1.replace('-', 'n'), seq2.replace('-', 'n'))
    except Exception as e:
        print(f"Error in pylcs.lcs_sequence_length with seq1: {seq1}, seq2: {seq2}")
        raise
    metrics['lcs_length'] = lcs
    metrics['normalized_lcs'] = lcs / min(len(seq1), len(seq2))
    # Normalize against ref sequence length
    metrics['normalized_lcs_ref'] = lcs / len(seq2)
    if timers is not None:
        timers['lcs_length'] += time.time() - start
    
    # PairwiseAligner score with default scoring
    if timers is not None:
        start = time.time()
    try:
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 5
        aligner.mismatch_score = -4
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        score = aligner.score(seq1, seq2)
    except Exception as e:
        print(f"Error in PairwiseAligner().score with seq1: {seq1}, seq2: {seq2}")
        raise
    metrics['pairwise_aligner'] = score
    # Normalized by the maximum possible score (all matches)
    perfect_score = aligner.match_score * min(len(seq1), len(seq2))
    metrics['normalized_pairwise_aligner'] = score / perfect_score if perfect_score != 0 else 0
    # Normalize against ref sequence length
    perfect_score_ref = aligner.match_score * len(seq2)
    metrics['normalized_pairwise_aligner_ref'] = score / perfect_score_ref if perfect_score_ref != 0 else 0
    
    # Perform global alignment (will be used for other metrics below)
    try:
        alignments = aligner.align(seq1, seq2)
        alignment = alignments[0]  # Get the first alignment
        aligned_seq1, aligned_seq2 = str(alignment[0]), str(alignment[1])
    except Exception as e:
        print(f"Error in PairwiseAligner().align with seq1: {seq1}, seq2: {seq2}")
        raise
    
    """ 
    # PairwiseAligner with BLASTN scoring
    if timers is not None:
        start_blastn = time.time()
    
    # Create a new aligner with BLASTN scoring
    try:
        blastn_aligner = PairwiseAligner(scoring="blastn")
        blastn_aligner.mode = 'global'
        
        # Use DNA-compatible sequences
        blastn_score = blastn_aligner.score(seq1_dna, seq2_dna)
        metrics['pairwise_aligner_blastn'] = blastn_score
        
        # Perform alignment with BLASTN scoring to get aligned sequences
        blastn_alignments = blastn_aligner.align(seq1_dna, seq2_dna)
        blastn_alignment = blastn_alignments[0]
        aligned_seq1_dna = str(blastn_alignment[0])
        aligned_seq2_dna = str(blastn_alignment[1])
        
        # Normalized BLASTN scores
        # For normalization, we need to estimate the perfect score
        # A rough estimate based on BLASTN match score (typically +1 for matches)
        blastn_perfect_score = min(len(seq1_dna), len(seq2_dna))  # Assuming match score of 1
        metrics['normalized_pairwise_aligner_blastn'] = blastn_score / blastn_perfect_score if blastn_perfect_score != 0 else 0
        metrics['normalized_pairwise_aligner_blastn_ref'] = blastn_score / len(seq2_dna) if len(seq2_dna) != 0 else 0
    except Exception as e:
        print(f"Error in PairwiseAligner(scoring='blastn') with seq1_dna: {seq1_dna}, seq2_dna: {seq2_dna}")
        raise
    
    if timers is not None:
        # Add BLASTN timing to pairwise aligner timing
        blastn_time = time.time() - start_blastn
        timers['pairwise_aligner'] += blastn_time
    
    # Calculate distances using Bio.Phylo.TreeConstruction.DistanceCalculator
    if timers is not None:
        start = time.time()


    # Use the already DNA-compatible aligned sequences from BLASTN alignment
    try:
        seq_record1 = SeqRecord(Seq(aligned_seq1_dna), id="seq1")
        seq_record2 = SeqRecord(Seq(aligned_seq2_dna), id="seq2")
        
        # Create a MultipleSeqAlignment object
        alignment_obj = MultipleSeqAlignment([seq_record1, seq_record2])
        
        # Calculate distances with different models
        for model in ['identity', 'blastn', 'trans']:
            calculator = DistanceCalculator(model)
            # Use the alignment object
            distance_matrix = calculator.get_distance(alignment_obj)
            # Extract the distance value
            distance_value = distance_matrix[1][0]  # The matrix is symmetric
            
            # Store the raw distance
            metrics[f'distance_{model}'] = distance_value
    except Exception as e:
        print(f"Error in DistanceCalculator with aligned_seq1_dna: {aligned_seq1_dna}, aligned_seq2_dna: {aligned_seq2_dna}")
        raise
    
    if timers is not None:
        # Add tree construction timing to pairwise aligner since they're related
        timers['pairwise_aligner'] += time.time() - start
    """
        
    return metrics

def process_record(input_record, ref_record, source_records, with_dtw):
    """Process a single record, computing metrics against both reference and source."""
    # Thread-local timers
    local_timers = {
        'levenshtein': 0,
        'dtw': 0,
        'jaro_winkler': 0,
        'lcs_length': 0,
        'pairwise_aligner': 0
    } if DEBUG_TIMING else None
    
    # Initialize variables that might be used in error handling
    ref_id = None
    source_record = None
    
    try:
        # Compare against reference alignment
        ref_metrics = compute_metrics(str(input_record.seq), str(ref_record.seq), with_dtw, local_timers)
        ref_metrics['id'] = input_record.id
        
        # Extract source reference ID
        source_metrics_entry = None
        ref_id = extract_reference_id(input_record.description)
        
        if ref_id:
            # Find the source record
            for source_id, record in source_records.items():
                if source_id.startswith(f"{ref_id} ") or source_id == ref_id:
                    source_record = record
                    break
            
            if source_record:
                # Compare against source data
                source_metrics_entry = compute_metrics(str(input_record.seq), str(source_record.seq), with_dtw, local_timers)
                source_metrics_entry['id'] = input_record.id
                source_metrics_entry['source_id'] = source_record.id
        else:
            print(f"Warning: No reference ID found for record {input_record.id}")
        
        return ref_metrics, source_metrics_entry, local_timers if DEBUG_TIMING else None
    except ValueError as e:
        if "sequence contains letters not in the alphabet" in str(e):
            # Print detailed information about the problematic sequences
            print(f"\nError processing record {input_record.id}:")
            print(f"Input sequence: {str(input_record.seq)}")
            print(f"Reference sequence: {str(ref_record.seq)}")
            if ref_id and source_record:
                print(f"Source sequence: {str(source_record.seq)}")
            print(f"Error message: {str(e)}")
            # Re-raise the exception to maintain the original error flow
            raise
        else:
            # Re-raise other ValueError exceptions
            raise

def main():
    args = parse_args()
    
    # Store all metrics
    reference_metrics = []
    source_metrics = []
    
    # Read input and reference alignment files
    input_records = list(SeqIO.parse(args.input, "fasta"))
    ref_align_records = list(SeqIO.parse(args.reference_alignment, "fasta"))
    
    # Verify that reference and input have the same number of records
    if len(input_records) != len(ref_align_records):
        print(f"Error: Input file has {len(input_records)} records, but reference file has {len(ref_align_records)} records")
        exit(1)
    
    # Verify that the headers match with reference alignment
    for i, (input_record, ref_record) in enumerate(zip(input_records, ref_align_records)):
        if input_record.id != ref_record.id:
            print(f"Error: Header mismatch at record {i} - Input: {input_record.id}, Reference: {ref_record.id}")
            exit(1)
    
    # Read all source data records for quick lookup
    source_records = {}
    for record in SeqIO.parse(args.source_data, "fasta"):
        source_records[record.id] = record
    
    # Combined metrics and timers
    combined_timers = {
        'levenshtein': 0,
        'dtw': 0,
        'jaro_winkler': 0,
        'lcs_length': 0,
        'pairwise_aligner': 0
    } if DEBUG_TIMING else None
    
    # Create pairs of input and reference records
    record_pairs = list(zip(input_records, ref_align_records))

    #record_pairs = record_pairs[:10]
    
    # Process records in parallel
    print(f"Using {args.threads} threads for processing...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(process_record, input_record, ref_record, source_records, args.with_dtw) 
                   for input_record, ref_record in record_pairs]
        
        for future in futures:
            ref_metric, source_metric, local_timers = future.result()
            reference_metrics.append(ref_metric)
            
            if source_metric:
                source_metrics.append(source_metric)
            
            # Aggregate timers
            if DEBUG_TIMING and local_timers:
                for key in combined_timers:
                    combined_timers[key] += local_timers[key]
    
    # Calculate average of metrics
    avg_ref_metrics = {}
    avg_source_metrics = {}
    
    metric_keys = [k for k in reference_metrics[0].keys() if k != 'id']
    
    for key in metric_keys:
        avg_ref_metrics[key] = np.mean([m[key] for m in reference_metrics])
        
    if source_metrics:
        source_metric_keys = [k for k in source_metrics[0].keys() if k not in ('id', 'source_id')]
        for key in source_metric_keys:
            avg_source_metrics[key] = np.mean([m[key] for m in source_metrics])
    
    # Output the results
    # 1. Stats file (averages)
    with open(f"{args.output_base_name}_compare_stats.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Metric', 'Reference_Avg', 'Source_Avg'])
        for key in metric_keys:
            writer.writerow([key, avg_ref_metrics[key], avg_source_metrics.get(key, 'N/A')])
    
    # 2. Reference scores file
    with open(f"{args.output_base_name}_reference_scores.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID'] + metric_keys)
        for metrics in reference_metrics:
            writer.writerow([metrics['id']] + [metrics[key] for key in metric_keys])
    
    # 3. Source scores file
    if source_metrics:
        with open(f"{args.output_base_name}_source_scores.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            source_metric_keys = [k for k in source_metrics[0].keys() if k not in ('id', 'source_id')]
            writer.writerow(['ID', 'Source_ID'] + source_metric_keys)
            for metrics in source_metrics:
                writer.writerow([metrics['id'], metrics['source_id']] + [metrics[key] for key in source_metric_keys])
    
    # Print averages to console
    print("Average Metrics:")
    print("Reference Comparisons:")
    for key in metric_keys:
        print(f"  {key}: {avg_ref_metrics[key]}")
    
    if source_metrics:
        print("\nSource Comparisons:")
        source_metric_keys = [k for k in source_metrics[0].keys() if k not in ('id', 'source_id')]
        for key in source_metric_keys:
            print(f"  {key}: {avg_source_metrics[key]}")
    
    # Print timing information if debugging is enabled
    if DEBUG_TIMING:
        print("\nTiming Information:")
        for metric, elapsed in combined_timers.items():
            print(f"  {metric}: {elapsed:.4f} seconds")

if __name__ == "__main__":
    main()

