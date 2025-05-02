#!/bin/bash

# Calculează acoperirea pentru fiecare set de date de test

for size in 1k 10k 50k 100k 200k 500k 1m 2m; do
    echo "Se calculează acoperirea pentru setul ${size}..."
    
    python3 compute_coverage.py -i "${size}.fasta"
    
    echo "----------------------------------------"
done
