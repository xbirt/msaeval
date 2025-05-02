#!/bin/bash

# Script pentru extragerea speciilor țintă din baza de date SILVA cu păstrarea metadatelor
# Cerințe: SeqKit instalat și accesibil în PATH

INPUT_FASTA="SILVA_138_SSURef_tax_silva.fasta"

# Definim toate categoriile de reads pentru care generăm subseturi
declare -a READ_CATEGORIES=(1000 10000 50000 100000 500000 1000000)

for NUM_READS in "${READ_CATEGORIES[@]}"; do
    declare -A SPECIES_PATTERNS
    OUTPUT_PREFIX="selected_species_${NUM_READS}"

    # Selectăm speciile în funcție de categoria de reads
    if [ "$NUM_READS" -le 1000 ]; then
        SPECIES_PATTERNS=(
            ["Escherichia_coli"]='\bEscherichia\s+coli\b'
            ["Staphylococcus_aureus"]='\bStaphylococcus\s+aureus\b'
        )
    elif [ "$NUM_READS" -le 10000 ]; then
        SPECIES_PATTERNS=(
            ["Escherichia_coli"]='\bEscherichia\s+coli\b'
            ["Staphylococcus_aureus"]='\bStaphylococcus\s+aureus\b'
            ["Bacillus_subtilis"]='\bBacillus\s+subtilis\b'
        )
    elif [ "$NUM_READS" -le 50000 ]; then
        SPECIES_PATTERNS=(
            ["Escherichia_coli"]='\bEscherichia\s+coli\b'
            ["Bacillus_subtilis"]='\bBacillus\s+subtilis\b'
            ["Staphylococcus_aureus_complex"]='Staphylococcus\s+(aureus|epidermidis)\b'
            ["Streptococcus_pneumoniae"]='Streptococcus\s+pneumoniae\b'
        )
    elif [ "$NUM_READS" -le 100000 ]; then
        SPECIES_PATTERNS=(
            ["Escherichia_coli"]='\bEscherichia\s+coli\b'
            ["Bacillus_subtilis"]='\bBacillus\s+subtilis\b'
            ["Staphylococcus_aureus_complex"]='Staphylococcus\s+(aureus|epidermidis)\b'
            ["Streptococcus_pneumoniae"]='Streptococcus\s+pneumoniae\b'
            ["Neisseria_meningitidis"]='\bNeisseria\s+meningitidis\b'
            ["Lactobacillus_iners"]='\bLactobacillus\s+iners\b'
        )
    elif [ "$NUM_READS" -le 500000 ]; then
        SPECIES_PATTERNS=(
            # Specii cu variație intra-genomică
            ["Escherichia_coli"]='\bEscherichia\s+coli\b'
            ["Bacillus_subtilis"]='\bBacillus\s+subtilis\b'
            
            # Complexe de specii apropiate
            ["Streptococcus_pneumoniae_complex"]='Streptococcus\s+(pneumoniae|mitis|oralis)\b'
            ["Staphylococcus_aureus_complex"]='Staphylococcus\s+(aureus|epidermidis)\b'
            
            # Specii cu regiuni hipervariabile
            ["Neisseria_meningitidis"]='\bNeisseria\s+meningitidis\b'
            ["Lactobacillus_iners"]='\bLactobacillus\s+iners\b'
            
            # Specii necultivate și ambiguități taxonomice
            ["Uncultured_bacteria"]='\buncultured\s+bacterium\b|\benvironmental\s+samples\b'
        #    ["Enterococcus_faecalis_conflict"]='Enterococcus\s+faecalis\b.*(conflict|SILVA|MIMt)'
        )
    else
        SPECIES_PATTERNS=(
            # Specii cu variație intra-genomică
            ["Escherichia_coli"]='\bEscherichia\s+coli\b'
            ["Bacillus_subtilis"]='\bBacillus\s+subtilis\b'
            
            # Complexe de specii apropiate
            ["Streptococcus_complex"]='Streptococcus\s\b'
            ["Staphylococcus_complex"]='Staphylococcus\s\b'
            
            # Specii cu regiuni hipervariabile
            ["Neisseria_meningitidis"]='\bNeisseria\s\b'
            ["Lactobacillus_iners"]='\bLactobacillus\s\b'
            
            # Specii necultivate și ambiguități taxonomice
            ["Uncultured_bacteria"]='\buncultured\s+bacterium\b|\benvironmental\s+samples\b'
        #    ["Enterococcus_faecalis_conflict"]='Enterococcus\s+faecalis\b.*(conflict|SILVA|MIMt)'
        )
    fi

    # Procesare pentru fiecare categorie
    echo -e "\nSe procesează $NUM_READS citiri..."
    OUTPUT_FILE="${OUTPUT_PREFIX}.fasta"
    OUTPUT_FILTERED_FILE="${OUTPUT_PREFIX}_filtered.fasta"
    OUTPUT_FILTERED_FILE_LENGTH="${OUTPUT_PREFIX}_filtered_length.fasta"
    
    > "$OUTPUT_FILE"
    
    for species in "${!SPECIES_PATTERNS[@]}"; do
        pattern=${SPECIES_PATTERNS[$species]}
        echo "Se extrage $species..."
        seqkit grep -n -r -p "$pattern" "$INPUT_FASTA" >> "$OUTPUT_FILE"
    done

    # Filtrare înregistrări sub 1400 bp
    seqkit seq -g -m 1400 "$OUTPUT_FILE" > "$OUTPUT_FILTERED_FILE_LENGTH"
    # Regex invers pe secvențe cu caracterele valide pentru nucleotide
    seqkit grep --by-seq --invert-match -r --pattern '[^ACTUGactug]' "$OUTPUT_FILTERED_FILE_LENGTH" > "$OUTPUT_FILTERED_FILE"
    
    echo "Rezultate pentru $NUM_READS citiri:"
    seqkit stats "$OUTPUT_FILTERED_FILE"
done

echo -e "\nToate fișierele au fost generate:"
ls -lh selected_species_*_filtered.fasta

