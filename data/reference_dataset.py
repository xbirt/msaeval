#!/usr/bin/env python3

# Construiește un fișier FASTA cu secvențele de referință unice corespunzătoare citirilor.
# Pentru fiecare referință găsită în citiri, vom avea o singură înregistrare în fișierul de ieșire.

import argparse
from typing import Set
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_references(header: str) -> Set[str]:
    """
    Extrage ID-urile de referință dintr-un antet FASTA.
    
    Args:
        header: Șirul de caractere care conține antetul FASTA
        
    Returns:
        Un set de șiruri de caractere reprezentând ID-urile de referință
    """
    # Împărțim header-ul în părți separate prin spații
    parts = header.split()
    # Căutăm partea care începe cu 'reference='
    for part in parts:
        if part.startswith('reference='):
            # Extragem referințele (tot ce e după 'reference=') și le împărțim folosind virgula ca separator
            return set(part[10:].split(','))
    return set()

def process_reads_file(reads_file: str) -> Set[str]:
    """
    Procesează fișierul de citiri FASTA și returnează setul de referințe unice.
    
    Args:
        reads_file: Calea către fișierul de citiri FASTA
        
    Returns:
        Un set de șiruri de caractere reprezentând toate referințele unice
    """
    references = set()
    # Parcurgem toate înregistrările din fișierul FASTA
    for record in SeqIO.parse(reads_file, "fasta"):
        # Extragem referințele din header și le adăugăm la set
        references.update(extract_references(record.description))
    return references

def process_reference_file(reference_file: str, references: Set[str]) -> dict:
    """
    Procesează fișierul sursă FASTA și returnează un dicționar cu secvențele corespunzătoare referințelor.
    
    Args:
        source_file: Calea către fișierul sursă FASTA
        references: Setul de referințe căutate
        
    Returns:
        Un dicționar care mapază ID-urile de referință la obiectele SeqRecord corespunzătoare
    """
    sequences = {}
    # Parcurgem toate înregistrările din fișierul sursă
    for record in SeqIO.parse(reference_file, "fasta"):
        # Extragem ID-ul de referință din header (tot ce e după '>' până la primul spațiu)
        ref_id = record.description.split()[0]
        # Dacă ID-ul este în setul de referințe căutate, adăugăm secvența în dicționar
        if ref_id in references:
            sequences[ref_id] = record
    return sequences

def main():
    """
    Funcția principală care coordonează procesarea fișierelor FASTA.
    """
    # Configurăm parser-ul de argumente din linia de comandă
    parser = argparse.ArgumentParser(description='Procesează fișiere FASTA și extrage secvențele de referință.')
    parser.add_argument('--reads', required=True, help='Fișierul de citiri FASTA de intrare')
    parser.add_argument('--reference', required=True, help='Fișierul de referință FASTA')
    parser.add_argument('--output', required=True, help='Fișierul FASTA de ieșire')
    
    args = parser.parse_args()
    
    # Procesăm fișierul de citiri pentru a obține referințele unice
    references = process_reads_file(args.reads)
    print(f"Am găsit {len(references)} referințe unice în fișierul de citiri")
    
    # Procesăm fișierul sursă pentru a obține secvențele corespunzătoare
    sequences = process_reference_file(args.reference, references)
    
    # Verificăm dacă toate referințele au fost găsite în fișierul sursă
    missing_refs = references - set(sequences.keys())
    if missing_refs:
        print("Eroare: Următoarele referințe nu au fost găsite în fișierul sursă:")
        for ref in sorted(missing_refs):
            print(f"  - {ref}")
        return
    
    # Scriem secvențele găsite în fișierul de ieșire
    output_records = []
    for ref in sorted(references):
        output_records.append(sequences[ref])
    
    SeqIO.write(output_records, args.output, "fasta")
    print(f"Am scris cu succes {len(references)} secvențe în {args.output}")

if __name__ == '__main__':
    main() 