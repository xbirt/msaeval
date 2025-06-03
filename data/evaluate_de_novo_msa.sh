#!/bin/bash

# Evaluate intrinsic scores
python score.py -i kalign.fasta --threads 8 --output-base-name kalign
python score.py -i pasta.fasta --threads 8 --output-base-name pasta
python score.py -i muscle.fasta --threads 8 --output-base-name muscle
python score.py -i clustalo.fasta --threads 8 --output-base-name clustalo
python score.py -i mafft.fasta --threads 8 --output-base-name mafft

# Evaluate scores against reference alignment and source data
python evaluate.py -i kalign.fasta --reference-alignment 10k-reference-full-alignment.fasta --source-data selected_species_10000_subset.fasta --threads 8 --dtw --output-base-name kalign
python evaluate.py -i pasta.fasta --reference-alignment 10k-reference-full-alignment.fasta --source-data selected_species_10000_subset.fasta --threads 8 --dtw --output-base-name pasta
python evaluate.py -i muscle.fasta --reference-alignment 10k-reference-full-alignment.fasta --source-data selected_species_10000_subset.fasta --threads 8 --dtw --output-base-name muscle
python evaluate.py -i clustalo.fasta --reference-alignment 10k-reference-full-alignment.fasta --source-data selected_species_10000_subset.fasta --threads 8 --dtw --output-base-name clustalo
python evaluate.py -i mafft.fasta --reference-alignment 10k-reference-full-alignment.fasta --source-data selected_species_10000_subset.fasta --threads 8 --dtw --output-base-name mafft

