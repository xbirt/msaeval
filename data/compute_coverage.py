#! /usr/bin/env python3

"""
Acest script calculează acoperirea (coverage) pentru fiecare referință dintr-un fișier FASTA.
Pentru fiecare secvență din fișier, extrage prima referință din header (după 'reference=') 
și numără de câte ori apare fiecare referință. La final, afișează numărul de apariții 
pentru fiecare referință și media acoperirii.

Utilizare:
    python compute_coverage.py --input <cale_catre_fisier_fasta> [--verbose]
    python compute_coverage.py -i <cale_catre_fisier_fasta> [-v]
"""

from Bio import SeqIO
from collections import defaultdict
import re
import argparse
import sys
import os

def get_first_reference(header):
    """Extrage prima referință din header-ul secvenței."""
    ref_match = re.search(r'reference=([^,\s]+)', header)
    return ref_match.group(1) if ref_match else None

def compute_coverage(fasta_file):
    """Calculează acoperirea pentru fiecare referință din fișierul FASTA."""
    coverage = defaultdict(int)
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        ref = get_first_reference(record.description)
        if ref:
            coverage[ref] += 1
    
    return coverage

def main():
    # Configurare parser pentru argumentele din linia de comandă
    parser = argparse.ArgumentParser(description='Calculează acoperirea pentru referințele dintr-un fișier FASTA')
    parser.add_argument('-i', '--input', '--reads', required=True, 
                       help='Calea către fișierul FASTA cu secvențele')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Afișează acoperirea detaliată pentru fiecare referință')
    
    try:
        args = parser.parse_args()
    except:
        print("Eroare: Parametrul de input este obligatoriu.", file=sys.stderr)
        print("Utilizare: python compute_coverage.py --input <cale_catre_fisier_fasta> [--verbose]", file=sys.stderr)
        sys.exit(1)
    
    # Verifică dacă fișierul există
    if not os.path.exists(args.input):
        print(f"Eroare: Fișierul '{args.input}' nu există.", file=sys.stderr)
        sys.exit(1)
    
    # Calculează acoperirea
    coverage = compute_coverage(args.input)
    
    # Afișează acoperirea detaliată dacă este specificat --verbose
    if args.verbose:
        print("Acoperire per referință:")
        for ref, cov in coverage.items():
            print(f"{ref}: {cov}")
    
    # Calculează și afișează media acoperirii
    if coverage:
        avg_coverage = sum(coverage.values()) / len(coverage)
        print(f"Media acoperirii: {avg_coverage:.2f}")
    else:
        print("Avertisment: Nu s-au găsit referințe în fișier.", file=sys.stderr)

if __name__ == "__main__":
    main()
