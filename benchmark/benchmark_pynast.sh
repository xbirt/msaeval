#!/bin/bash

# Download the SILVA 138 SSU Ref aligned databases if needed
if [ ! -f "SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta" ]; then
    if [ ! -f "SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz" ]; then
        wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz
    fi
    gunzip SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz
fi

if [ ! -f "SILVA_138_SSURef_tax_silva_full_align_trunc.fasta" ]; then
    if [ ! -f "SILVA_138_SSURef_tax_silva_full_align_trunc.fasta.gz" ]; then
        wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_tax_silva_full_align_trunc.fasta.gz
    fi
    gunzip SILVA_138_SSURef_tax_silva_full_align_trunc.fasta.gz
fi

./benchmark.sh "pynast" "-i 1k.fasta -t SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 10k.fasta -t SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 50k.fasta -t SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 100k.fasta -t SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 500k.fasta -t SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 1m.fasta -t SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 2m.fasta -t SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta"

./benchmark.sh "pynast" "-i 1k.fasta -t SILVA_138_SSURef_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 10k.fasta -t SILVA_138_SSURef_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 50k.fasta -t SILVA_138_SSURef_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 100k.fasta -t SILVA_138_SSURef_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 500k.fasta -t SILVA_138_SSURef_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 1m.fasta -t SILVA_138_SSURef_tax_silva_full_align_trunc.fasta"
./benchmark.sh "pynast" "-i 2m.fasta -t SILVA_138_SSURef_tax_silva_full_align_trunc.fasta"