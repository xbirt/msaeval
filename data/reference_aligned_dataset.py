#!/usr/bin/env python3

# Crează un fișier FASTA cu secvențele aliniate la referință, cu o corespondență 1:1 între citiri și referințe.
# Pentru fiecare citire, vom avea câte o secvență din fișierul de referință în fișierul de ieșire, ajustată la poziția ampliconului.
# Chimerele sunt de asemenea procesate în acest script.

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys
import argparse
import os
from typing import List, Tuple, Dict

def parse_header(header: str) -> Tuple[str, List[str], List[Tuple[int, int]], Tuple[int, int]]:
    """Parsează header-ul FASTA pentru a extrage ID-ul, referințele, pozițiile ampliconului și poziția secvenței.
    
    Args:
        header: String-ul header-ului FASTA
        
    Returns:
        Tuple care conține:
        - ID-ul secvenței
        - Lista de ID-uri de referință
        - Lista de tuple (start, end) pentru fiecare referință
        - Tuple (start, end) pentru poziția finală
    """
    # Extragem ID-ul
    id_match = re.match(r'(\S+)', header)
    if not id_match:
        raise ValueError(f"Format invalid al header-ului: {header}")
    seq_id = id_match.group(1)
    
    # Extragem referințele
    ref_match = re.search(r'reference=([^\s,]+)', header)
    if not ref_match:
        raise ValueError(f"Nu s-a găsit referință în header: {header}")
    references = ref_match.group(1).split(',')
    
    # Extragem pozițiile ampliconului
    amp_match = re.search(r'amplicon=([^\s,]+)', header)
    if not amp_match:
        raise ValueError(f"Nu s-au găsit poziții de amplicon în header: {header}")
    amplicons = []
    for pos in amp_match.group(1).split(','):
        start, end = map(int, pos.split('..'))
        amplicons.append((start, end))
    
    # Extragem poziția secvenței
    pos_match = re.search(r'position=(\d+)\.\.(\d+)', header)
    if not pos_match:
        raise ValueError(f"Nu s-a găsit poziție în header: {header}")
    position = (int(pos_match.group(1)), int(pos_match.group(2)))
    
    return seq_id, references, amplicons, position

def extract_sequence_from_reference(ref_id: str, start: int, end: int, 
                                  reference_records: Dict[str, SeqRecord]) -> str:
    """Extrage o subsecvență dintr-o secvență de referință.
    
    Args:
        ref_id: ID-ul referinței de căutat
        start: Poziția de început (bazată pe 1)
        end: Poziția de sfârșit (bazată pe 1)
        reference_records: Dicționar cu secvențele de referință
        
    Returns:
        Secvența extrasă ca string
    """
    # Căutăm înregistrarea de referință care începe cu ID-ul dat
    for ref_record in reference_records.values():
        if ref_record.id.startswith(ref_id):
            # Convertim la indexare bazată pe 0 pentru Python
            return str(ref_record.seq[start-1:end])
    raise ValueError(f"Referința {ref_id} nu a fost găsită în fișierul de referință")

def process_fasta(reads_file: str, reference_file: str, output_file: str):
    """Procesează fișierele FASTA de intrare și generează output-ul aliniat la referință.
    
    Args:
        reads_file: Calea către fișierul FASTA cu citiri
        reference_file: Calea către fișierul FASTA de referință
        output_file: Calea către fișierul FASTA de ieșire
    """
    # Citim secvențele de referință
    reference_records = {record.id: record for record in SeqIO.parse(reference_file, "fasta")}
    
    # Procesăm citirile și generăm output-ul
    output_records = []
    for record in SeqIO.parse(reads_file, "fasta"):
        try:
            # Parsează header-ul
            seq_id, references, amplicons, position = parse_header(record.description)
            
            # Extragem secvențele din referințe
            sequences = []
            for ref_id, (start, end) in zip(references, amplicons):
                seq = extract_sequence_from_reference(ref_id, start, end, reference_records)
                sequences.append(seq)
            
            # Concatenăm secvențele și extragem poziția finală
            full_sequence = ''.join(sequences)
            final_sequence = full_sequence[position[0]-1:position[1]]
            
            # Creăm înregistrarea de output
            output_record = SeqRecord(
                Seq(final_sequence),
                id=seq_id,
                description=record.description[1:]  # Eliminăm '>' din descrierea originală
            )
            output_records.append(output_record)
            
        except ValueError as e:
            print(f"Eroare la procesarea înregistrării {record.id}: {str(e)}", file=sys.stderr)
    
    # Scriem output-ul
    SeqIO.write(output_records, output_file, "fasta")

if __name__ == "__main__":
    # Configurăm parser-ul de argumente din linia de comandă
    parser = argparse.ArgumentParser(description='Procesează fișiere FASTA și extrage secvențele de referință.')
    parser.add_argument('--reads', required=True, help='Fișierul de citiri FASTA de intrare')
    parser.add_argument('--reference', required=True, help='Fișierul de referință FASTA')
    parser.add_argument('--output', required=True, help='Fișierul FASTA de ieșire')
    
    args = parser.parse_args()
    
    # Verificăm dacă fișierele de intrare există
    if not os.path.exists(args.reads):
        print(f"Eroare: Fișierul de citiri '{args.reads}' nu există.", file=sys.stderr)
        sys.exit(1)
        
    if not os.path.exists(args.reference):
        print(f"Eroare: Fișierul de referință '{args.reference}' nu există.", file=sys.stderr)
        sys.exit(1)
    
    process_fasta(args.reads, args.reference, args.output)