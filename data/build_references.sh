#!/bin/bash

for size in 1k 10k 50k 100k 200k 500k 1m 2m; do
    subset_size=$(case $size in
        "1k") echo "1000";;
        "10k") echo "10000";;
        "50k") echo "50000";;
        "100k") echo "100000";;
        "200k") echo "100000";;
        "500k") echo "500000";;
        "1m") echo "1000000";;
        "2m") echo "1000000";;
    esac)

    # Construim referințele pentru toate datele de test - doar referințe unice
    if [ ! -f "${size}-reference.fasta" ] || [ ! -s "${size}-reference.fasta" ]; then
        echo "Se construiește setul de referințe unice pentru ${size}..."
        python3 reference_dataset.py \
            --reads ${size}.fasta \
            --reference selected_species_${subset_size}_subset.fasta \
            --output ${size}-reference.fasta
    fi

    # Construim referințele pentru toate datele de test - câte o referință pentru fiecare citire
    if [ ! -f "${size}-reference-full.fasta" ] || [ ! -s "${size}-reference-full.fasta" ]; then
        echo "Se construiește setul complet de referințe pentru ${size}..."
        python3 reference_aligned_dataset.py \
            --reads ${size}.fasta \
            --reference selected_species_${subset_size}_subset.fasta \
            --output ${size}-reference-full.fasta
    fi
done
