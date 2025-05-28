#!/usr/bin/env python3

import os
import csv
import sys
from collections import OrderedDict
import jellyfish
import Levenshtein
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from fastdtw import fastdtw

"""
This script aggregates and compares alignment and consensus metrics for multiple sequence alignment (MSA) tools.

For each tool (name in NAMES), it:
  1. Reads stats and compare_stats CSVs, storing all metrics.
  2. For each consensus sequence variant (50/70, with/without gaps):
     - Compares the tool's consensus to the reference consensus using:
         * Jaro-Winkler similarity
         * Levenshtein distance (raw and normalized)
         * (For no_gaps only) Global pairwise alignment score (Biopython PairwiseAligner, NUC44-like scoring) and normalized score
  3. Checks that all metrics are present for all tools.
  4. Outputs a centralized CSV with all metrics as rows and all tools as columns.

The output file is 'centralized_de_novo_msa_alignment_metrics.csv'.
"""

NAMES = [
    '10k-reference-full-alignment',
    'mafft1',
    'clustalo1',
    'kalign1',
    'pasta1',
    'muscle1',
]

CONSENSUS_FILES = [
    '_consensus_50_no_gaps.fasta',
    '_consensus_50_with_gaps.fasta',
    '_consensus_70_no_gaps.fasta',
    '_consensus_70_with_gaps.fasta',
]

REFERENCE_CONSENSUS = [
    '10k-reference-full-alignment_consensus_50_no_gaps.fasta',
    '10k-reference-full-alignment_consensus_50_with_gaps.fasta',
    '10k-reference-full-alignment_consensus_70_no_gaps.fasta',
    '10k-reference-full-alignment_consensus_70_with_gaps.fasta',
]

STATS_SUFFIX = '_stats.csv'
COMPARE_STATS_SUFFIX = '_compare_stats.csv'

# Step 1: Read stats files
def read_stats_file(filename):
    stats = OrderedDict()
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if ',' in line:
                key, value = line.split(',', 1)
                stats[key.strip()] = value.strip()
    return stats

# Step 2: Read compare_stats files
def read_compare_stats_file(filename):
    stats = OrderedDict()
    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if len(row) < 3:
                continue
            metric, ref, src = row[0].strip(), row[1].strip(), row[2].strip()
            stats[f'{metric} vs ref'] = ref
            stats[f'{metric} vs src'] = src
    return stats

# Step 3: Read consensus sequence from fasta
def read_single_fasta_sequence(filename):
    try:
        with open(filename, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                return str(record.seq)
    except Exception as e:
        print(f'Error reading {filename}: {e}', file=sys.stderr)
        sys.exit(1)
    print(f'Error: No sequence found in {filename}', file=sys.stderr)
    sys.exit(1)

# Helper to format metric names
def format_consensus_metric_name(cons_suffix, metric_name):
    # cons_suffix: '_consensus_50_no_gaps.fasta' -> 'Consensus 50 no gaps [metric_name]'
    base = cons_suffix.replace('_consensus_', '').replace('.fasta', '')
    parts = base.split('_')
    threshold = parts[0]
    gaps = 'with gaps' if 'with' in parts else 'no gaps'
    return f'Consensus {threshold} {gaps} {metric_name}'

# Remove ambiguous nucleotides from a sequence
def remove_ambiguous(seq):
    return ''.join([c for c in seq if c not in ('n', 'N')])

# Perform global pairwise alignment and return the score and normalized score
def global_alignment_score(seq1, seq2):
    aligner = PairwiseAligner()
    # Use NUC44 scoring matrix values: match=5, mismatch=-4, open=-10, extend=-0.5
    aligner.mode = 'global'
    match = 5
    mismatch = -4
    gap_open = -10
    gap_extend = -0.5
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    score = aligner.score(seq1, seq2)
    max_len = max(len(seq1), len(seq2))
    max_score = match * max_len
    min_score = min(gap_open, mismatch) * max_len
    norm = (score - min_score) / (max_score - min_score) if max_score != min_score else 1.0
    return score, norm

# DTW distance, average per-step, and normalized similarity
# Use fastdtw, which returns the path
# Normalization 1: average per-step cost (dtw / path_length)
# Normalization 2: similarity = 1 - (dtw / max_possible_cost), where max_possible_cost = path_length * max_per_step_cost
# For nucleotide ordinals, max_per_step_cost = max(abs(ord(a)-ord(b))) for a,b in 'ACGTN-'
def dtw_distance(seq1, seq2):
    s1 = [ord(c) for c in seq1.upper()]
    s2 = [ord(c) for c in seq2.upper()]
    dist, path = fastdtw(s1, s2)
    path_length = len(path)
    avg_per_step = dist / path_length if path_length > 0 else float('nan')
    # Compute max possible per-step cost for nucleotides
    alphabet = 'ACGTN-'
    max_per_step_cost = max(abs(ord(a) - ord(b)) for a in alphabet for b in alphabet)
    max_possible_cost = path_length * max_per_step_cost
    norm_similarity = 1.0 - (dist / max_possible_cost) if max_possible_cost > 0 else float('nan')
    return dist, avg_per_step, norm_similarity

def main():
    all_stats = {}
    metric_order = []
    metric_set = set()
    consensus_jw = {name: OrderedDict() for name in NAMES}
    consensus_lev = {name: OrderedDict() for name in NAMES}
    consensus_levnorm = {name: OrderedDict() for name in NAMES}
    consensus_align = {name: OrderedDict() for name in NAMES}
    consensus_alignnorm = {name: OrderedDict() for name in NAMES}
    #consensus_dtw = {name: OrderedDict() for name in NAMES}
    #consensus_dtw_avg = {name: OrderedDict() for name in NAMES}
    #consensus_dtw_norm = {name: OrderedDict() for name in NAMES}
    # Read all files and collect metrics
    for idx, name in enumerate(NAMES):
        stats_file = f'{name}{STATS_SUFFIX}'
        compare_stats_file = f'{name}{COMPARE_STATS_SUFFIX}'
        if not os.path.exists(stats_file):
            print(f'Error: File not found: {stats_file}', file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(compare_stats_file):
            print(f'Error: File not found: {compare_stats_file}', file=sys.stderr)
            sys.exit(1)
        stats = read_stats_file(stats_file)
        compare_stats = read_compare_stats_file(compare_stats_file)
        combined = OrderedDict()
        # Maintain order: stats first, then compare_stats
        for k, v in stats.items():
            combined[k] = v
            if idx == 0:
                metric_order.append(k)
                metric_set.add(k)
        for k, v in compare_stats.items():
            combined[k] = v
            if idx == 0:
                metric_order.append(k)
                metric_set.add(k)
        all_stats[name] = combined
    # Consensus Jaro-Winkler, Levenshtein, alignment score for no_gaps, and DTW
    consensus_jw_metric_names = []
    consensus_lev_metric_names = []
    consensus_levnorm_metric_names = []
    consensus_align_metric_names = []
    consensus_alignnorm_metric_names = []
    #consensus_dtw_metric_names = []
    #consensus_dtw_avg_metric_names = []
    #consensus_dtw_norm_metric_names = []
    for i, cons_suffix in enumerate(CONSENSUS_FILES):
        ref_file = REFERENCE_CONSENSUS[i]
        if not os.path.exists(ref_file):
            print(f'Error: Reference consensus file not found: {ref_file}', file=sys.stderr)
            sys.exit(1)
        ref_seq = read_single_fasta_sequence(ref_file)
        jw_metric_name = format_consensus_metric_name(cons_suffix, 'Jaro-Winkler similarity')
        lev_metric_name = format_consensus_metric_name(cons_suffix, 'Levenshtein distance')
        levnorm_metric_name = format_consensus_metric_name(cons_suffix, 'normalized Levenshtein similarity')
        #dtw_metric_name = format_consensus_metric_name(cons_suffix, 'DTW distance')
        #dtw_avg_metric_name = format_consensus_metric_name(cons_suffix, 'average per-step DTW')
        #dtw_norm_metric_name = format_consensus_metric_name(cons_suffix, 'normalized DTW similarity')
        consensus_jw_metric_names.append(jw_metric_name)
        consensus_lev_metric_names.append(lev_metric_name)
        consensus_levnorm_metric_names.append(levnorm_metric_name)
        #consensus_dtw_metric_names.append(dtw_metric_name)
        #consensus_dtw_avg_metric_names.append(dtw_avg_metric_name)
        #consensus_dtw_norm_metric_names.append(dtw_norm_metric_name)
        # Only for no_gaps
        if 'no_gaps' in cons_suffix:
            align_metric_name = format_consensus_metric_name(cons_suffix, 'pairwise alignment score')
            alignnorm_metric_name = format_consensus_metric_name(cons_suffix, 'normalized pairwise alignment score')
            consensus_align_metric_names.append(align_metric_name)
            consensus_alignnorm_metric_names.append(alignnorm_metric_name)
        for name in NAMES:
            cons_file = f'{name}{cons_suffix}'
            if not os.path.exists(cons_file):
                print(f'Error: Consensus file not found: {cons_file}', file=sys.stderr)
                sys.exit(1)
            seq = read_single_fasta_sequence(cons_file)
            # Jaro-Winkler similarity
            jw = jellyfish.jaro_winkler_similarity(seq, ref_seq)
            consensus_jw[name][jw_metric_name] = jw
            # Levenshtein distance
            lev = Levenshtein.distance(seq, ref_seq)
            consensus_lev[name][lev_metric_name] = lev
            # Normalized Levenshtein similarity (reuse lev)
            if len(seq) == 0 and len(ref_seq) == 0:
                levnorm = 1.0
            else:
                max_len = max(len(seq), len(ref_seq))
                levnorm = 1.0 - (lev / max_len)
            consensus_levnorm[name][levnorm_metric_name] = levnorm
            # DTW distance, average per-step, and normalized similarity
            #dtw_val, dtw_avg, dtw_norm = dtw_distance(seq, ref_seq)
            #consensus_dtw[name][dtw_metric_name] = dtw_val
            #consensus_dtw_avg[name][dtw_avg_metric_name] = dtw_avg
            #consensus_dtw_norm[name][dtw_norm_metric_name] = dtw_norm
            # Pairwise alignment score for no_gaps
            if 'no_gaps' in cons_suffix:
                seq_clean = remove_ambiguous(seq)
                ref_seq_clean = remove_ambiguous(ref_seq)
                align_score, align_norm = global_alignment_score(seq_clean, ref_seq_clean)
                consensus_align[name][align_metric_name] = align_score
                consensus_alignnorm[name][alignnorm_metric_name] = align_norm
    # Ensure all metrics are present in all structures
    for name in NAMES:
        missing = set(metric_order) - set(all_stats[name].keys())
        if missing:
            print(f'Error: {name} is missing metrics: {missing}', file=sys.stderr)
            sys.exit(1)
    # Write output
    with open('centralized_de_novo_msa_alignment_metrics.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['metric'] + NAMES)
        for metric in metric_order:
            row = [metric] + [all_stats[name][metric] for name in NAMES]
            writer.writerow(row)
        # Write consensus Jaro-Winkler similarity
        for metric in consensus_jw_metric_names:
            row = [metric] + [consensus_jw[name][metric] for name in NAMES]
            writer.writerow(row)
        # Write consensus Levenshtein distance
        for metric in consensus_lev_metric_names:
            row = [metric] + [consensus_lev[name][metric] for name in NAMES]
            writer.writerow(row)
        # Write consensus normalized Levenshtein similarity
        for metric in consensus_levnorm_metric_names:
            row = [metric] + [consensus_levnorm[name][metric] for name in NAMES]
            writer.writerow(row)
        # Write consensus DTW distance
        #for metric in consensus_dtw_metric_names:
        #    row = [metric] + [consensus_dtw[name][metric] for name in NAMES]
        #    writer.writerow(row)
        # Write consensus average per-step DTW
        #for metric in consensus_dtw_avg_metric_names:
        #    row = [metric] + [consensus_dtw_avg[name][metric] for name in NAMES]
        #    writer.writerow(row)
        # Write consensus normalized DTW similarity
        #for metric in consensus_dtw_norm_metric_names:
        #    row = [metric] + [consensus_dtw_norm[name][metric] for name in NAMES]
        #    writer.writerow(row)
        # Write consensus pairwise alignment scores for no_gaps
        for metric in consensus_align_metric_names:
            row = [metric] + [consensus_align[name].get(metric, '') for name in NAMES]
            writer.writerow(row)
        # Write consensus normalized pairwise alignment scores for no_gaps
        for metric in consensus_alignnorm_metric_names:
            row = [metric] + [consensus_alignnorm[name].get(metric, '') for name in NAMES]
            writer.writerow(row)

if __name__ == '__main__':
    main()
