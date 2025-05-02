#!/bin/bash

# Acest script procesează fișierele output curesim cu ajutorul porechop pentru eliminarea adaptoarelor și a chimerelor.
# De asemenea, filtrează rezultatele obținute, păstrând doar citirile de peste 1000 bp.

# Procesăm fiecare fișier curesim-*-reads.fastq
for input_file in curesim-*-reads.fastq; do
    # Extragem numele de bază (e.g., "1k" din "curesim-1k-reads.fastq")
    prefix=$(basename "$input_file" | cut -d- -f2)
    
    echo "Se procesează $prefix..."
    
    # Pas 1: Executăm Porechop cu procesare de chimere (split)
    porechop -i "$input_file" \
        -o "${prefix}-porechop.fastq" \
        --format fastq \
        --threads 20 \
        > "${prefix}-porechop.log" 2>&1
    
    # Pas 2: Eliminăm citirile <1000 bp
    seqkit seq -m 1000 -j 20 "${prefix}-porechop.fastq" -o "${prefix}-porechop-filtered.fastq"

    # Pas 3: Convertim fișierul rezultat FASTQ în FASTA și le păstrăm pe ambele
    seqkit fq2fa "${prefix}-porechop-filtered.fastq" -o "${prefix}-porechop-filtered.fasta" -j 20
    
    # Pas 4: Ștergem fișierul porechop care nu mai este necesar
    #rm "${prefix}-porechop.fastq"
done

echo "Procesarea a fost finalizată."