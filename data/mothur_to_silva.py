#! /usr/bin/env python3

from Bio import SeqIO
import os
from collections import defaultdict
from Bio import Align

# Identify the accession mapping from the Mothur-curated databases to the Silva database

# File name constants
MOTHUR_ALIGNMENT = "silva.nr_v138_2.align"
#MOTHUR_ALIGNMENT = "silva.seed_v138_2.align"
SILVA_DB = "SILVA_138.2_SSURef_tax_silva.fasta"
OUTPUT_FILE = "mothur2silva-nr99.txt"
ORPHANS_FILE = "orphans-nr99.txt"
#OUTPUT_FILE = "mothur2silva-seed.txt"
#ORPHANS_FILE = "orphans-seed.txt"

def extract_reference(header):
    """Extract reference from header (string before first space)"""
    return header.split()[0]

def extract_prefix(reference):
    """Extract prefix from reference (part before first dot)"""
    return reference.split('.')[0]

def clean_sequence(sequence):
    """Remove dots and dashes from sequence"""
    return str(sequence).replace('.', '').replace('-', '')

def find_best_alignment_match(mothur_group, silva_group):
    """Find the best alignment match between mothur and silva sequences"""
    aligner = Align.PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    best_score = -float('inf')
    best_mothur_ref = None
    best_silva_ref = None
    
    for mothur_ref, mothur_seq in mothur_group.items():
        mothur_seq_rna = mothur_seq.replace('T', 'U')
        
        for silva_ref, silva_seq in silva_group.items():
            # Perform pairwise alignment
            alignments = aligner.align(mothur_seq_rna, silva_seq)
            if alignments:
                score = alignments[0].score
                if score > best_score:
                    best_score = score
                    best_mothur_ref = mothur_ref
                    best_silva_ref = silva_ref
    
    return best_mothur_ref, best_silva_ref, best_score

def find_best_alignment_match_for_sequence(mothur_seq, silva_group):
    """Find the best alignment match for a single mothur sequence against all silva sequences in the group"""
    aligner = Align.PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    best_score = -float('inf')
    best_silva_ref = None
    
    mothur_seq_rna = mothur_seq.replace('T', 'U')
    
    for silva_ref, silva_seq in silva_group.items():
        try:
            # Perform pairwise alignment
            alignments = aligner.align(mothur_seq_rna, silva_seq)
            
            # Get the score without checking if alignments exist (to avoid OverflowError)
            score = alignments.score
            if score > best_score:
                best_score = score
                best_silva_ref = silva_ref
                
        except Exception as e:
            # Handle any other alignment errors
            print(f"Warning: Alignment error for {silva_ref}: {e}")
            print('mothur_seq_rna: ', mothur_seq_rna)
            print('silva_seq: ', silva_seq)
            raise e
    
    return best_silva_ref, best_score

def main():
    # Step 1: Read MOTHUR_ALIGNMENT and build mothur dictionaries
    print(f"Reading {MOTHUR_ALIGNMENT}...")
    mothur_by_prefix = defaultdict(dict)  # prefix -> {reference: sequence}
    mothur_sequences = {}  # reference -> sequence (for orphan tracking)
    
    try:
        with open(MOTHUR_ALIGNMENT, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                reference = extract_reference(record.id)
                prefix = extract_prefix(reference)
                cleaned_sequence = clean_sequence(record.seq)
                
                mothur_by_prefix[prefix][reference] = cleaned_sequence
                mothur_sequences[reference] = cleaned_sequence
        
        print(f"Loaded {len(mothur_sequences)} sequences from {MOTHUR_ALIGNMENT}")
        print(f"Organized into {len(mothur_by_prefix)} prefix groups")
    except FileNotFoundError:
        print(f"Error: {MOTHUR_ALIGNMENT} not found!")
        return
    
    # Step 2: Read SILVA_DB and build silva dictionaries
    print(f"Reading {SILVA_DB}...")
    silva_by_prefix = defaultdict(dict)  # prefix -> {reference: sequence}
    
    try:
        with open(SILVA_DB, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                reference = extract_reference(record.id)
                prefix = extract_prefix(reference)
                sequence = str(record.seq)
                
                silva_by_prefix[prefix][reference] = sequence
        
        total_silva_sequences = sum(len(refs) for refs in silva_by_prefix.values())
        print(f"Loaded {total_silva_sequences} sequences from {SILVA_DB}")
        print(f"Organized into {len(silva_by_prefix)} prefix groups")
    except FileNotFoundError:
        print(f"Error: {SILVA_DB} not found!")
        return
    
    # Step 3: Find matches by comparing sequences within matching prefix groups
    print("Finding matches using prefix-based optimization...")
    refdict = {}
    matched_mothur_refs = set()
    
    # Find common prefixes
    common_prefixes = set(mothur_by_prefix.keys()) & set(silva_by_prefix.keys())
    print(f"Found {len(common_prefixes)} common prefixes to search")
    
    forced_matches = 0
    alignment_matches = 0
    
    for prefix in common_prefixes:
        mothur_group = mothur_by_prefix[prefix]
        silva_group = silva_by_prefix[prefix]
        
        # Track which mothur sequences in this group have been matched
        group_matched_mothur_refs = set()
        
        # First pass: Try direct matches (exact and containment)
        for mothur_ref, mothur_seq in mothur_group.items():
            if mothur_ref in group_matched_mothur_refs:
                continue
                
            mothur_seq_rna = mothur_seq.replace('T', 'U')
            
            for silva_ref, silva_seq in silva_group.items():
                # Check if mothur sequence is contained in silva sequence
                if mothur_seq_rna in silva_seq:
                    refdict[mothur_ref] = silva_ref
                    matched_mothur_refs.add(mothur_ref)
                    group_matched_mothur_refs.add(mothur_ref)
                    break
        
        # Second pass: If no direct matches found, use alignment
        unmatched_mothur = {ref: seq for ref, seq in mothur_group.items() 
                           if ref not in group_matched_mothur_refs}
        
        if unmatched_mothur and silva_group:
            # If both groups have single entries, force the match
            if len(unmatched_mothur) == 1 and len(silva_group) == 1:
                mothur_ref = list(unmatched_mothur.keys())[0]
                silva_ref = list(silva_group.keys())[0]
                refdict[mothur_ref] = silva_ref
                matched_mothur_refs.add(mothur_ref)
                forced_matches += 1
            else:
                # Use pairwise alignment to find best match for each mothur sequence
                for mothur_ref, mothur_seq in unmatched_mothur.items():
                    best_silva_ref, best_score = find_best_alignment_match_for_sequence(
                        mothur_seq, silva_group)
                    
                    if best_silva_ref:
                        refdict[mothur_ref] = best_silva_ref
                        matched_mothur_refs.add(mothur_ref)
                        alignment_matches += 1
                
                if not any(ref in matched_mothur_refs for ref in unmatched_mothur.keys()):
                    print(f"No alignment matches found for prefix {prefix}: {list(unmatched_mothur.keys())}")
    
    print(f"Found {len(refdict)} matches between the two databases")
    print(f"Of which {forced_matches} were forced matches for single-entry prefix groups")
    print(f"And {alignment_matches} were found using pairwise alignment")
    
    # Step 4: Save mappings to OUTPUT_FILE
    print(f"Saving mappings to {OUTPUT_FILE}...")
    with open(OUTPUT_FILE, "w") as outfile:
        for mothur_ref, silva_ref in refdict.items():
            outfile.write(f"{mothur_ref}\t{silva_ref}\n")
    
    # Step 5: Save orphans to ORPHANS_FILE
    print(f"Saving orphans to {ORPHANS_FILE}...")
    orphan_refs = set(mothur_sequences.keys()) - matched_mothur_refs
    
    with open(ORPHANS_FILE, "w") as outfile:
        for orphan_ref in orphan_refs:
            outfile.write(f"{orphan_ref}\n")
    
    # Step 6: Save FASTA files for orphaned sequences by prefix
    print("Saving FASTA files for orphaned sequences by prefix...")
    
    # Group orphan references by prefix
    orphan_by_prefix = defaultdict(set)
    for orphan_ref in orphan_refs:
        prefix = extract_prefix(orphan_ref)
        orphan_by_prefix[prefix].add(orphan_ref)
    
    fasta_files_created = 0
    for prefix, orphan_refs_in_prefix in orphan_by_prefix.items():
        if orphan_refs_in_prefix:  # Only create files if there are orphans for this prefix
            # Create combined FASTA file for this prefix
            fasta_filename = f"orphans-{prefix}.fasta"
            with open(fasta_filename, "w") as fasta_file:
                # Write orphaned mothur sequences for this prefix
                if prefix in mothur_by_prefix:
                    for ref, seq in mothur_by_prefix[prefix].items():
                        if ref in orphan_refs_in_prefix:
                            fasta_file.write(f">{ref}\n{seq.replace('T', 'U')}\n")
                
                # Write all silva sequences for this prefix
                if prefix in silva_by_prefix:
                    for ref, seq in silva_by_prefix[prefix].items():
                        fasta_file.write(f">{ref}\n{seq}\n")
            
            fasta_files_created += 1
            #print(f"  Created {fasta_filename}")
    
    print(f"Orphaned references (no match in Silva): {len(orphan_refs)}")
    print(f"Created {fasta_files_created} FASTA files for {len(orphan_by_prefix)} orphaned prefixes")
    print("Processing complete!")
    print(f"Results saved to:")
    print(f"  - {OUTPUT_FILE} ({len(refdict)} mappings)")
    print(f"  - {ORPHANS_FILE} ({len(orphan_refs)} orphans)")
    print(f"  - {fasta_files_created} orphan FASTA files")

if __name__ == "__main__":
    main()

