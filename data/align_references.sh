#!/bin/bash

# Inițial proiectat să folosească mafft-xinsi cu extensia mxscarna, însă spațiul de disc necesar este prea mare.

# Construim aliniamente pentru toate referințele
for size in 1k 10k 50k 100k 200k 500k 1m 2m; do
    input_file="$size-reference.fasta"
    output_file="$size-reference-alignment.fasta"
    
    # Dacă fișierul de ieșire există și are dimensiune nenulă, sărim acest pas
    if [ -f "$output_file" ] && [ -s "$output_file" ]; then
        echo "Sărim peste ${output_file} - fișierul de aliniament există deja și are dimensiune nenulă."
        continue
    fi

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
    mafft --thread $(nproc) \
        $input_file \
        > $output_file
#    mafft-xinsi --thread $(nproc) \
#        --scarnapair \
#        $input_file \
#        > $output_file
done
