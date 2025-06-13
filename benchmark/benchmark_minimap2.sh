#!/bin/bash

THREADS=$(nproc)

# Create the minimap2 index for the SILVA database if needed
if [ ! -f "SILVA_138_SSURef_tax_silva.mmi" ]; then
    ./benchmark.sh "minimap2/minimap2" "-x map-ont -d SILVA_138_SSURef_tax_silva.mmi SILVA_138_SSURef_tax_silva.fasta"
fi

if [ ! -f "SILVA_138_SSURef_NR99_tax_silva.mmi" ]; then
    ./benchmark.sh "minimap2/minimap2" "-x map-ont -d SILVA_138_SSURef_NR99_tax_silva.mmi SILVA_138_SSURef_NR99_tax_silva.fasta"
fi

./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_tax_silva.mmi 1k.fasta > 1k-alignment.paf"
./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_tax_silva.mmi 10k.fasta > 10k-alignment.paf"
./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_tax_silva.mmi 50k.fasta > 50k-alignment.paf"
./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_tax_silva.mmi 100k.fasta > 100k-alignment.paf"
./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_tax_silva.mmi 500k.fasta > 500k-alignment.paf"
#./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_tax_silva.mmi 1m.fasta > 1m-alignment.paf"
#./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_tax_silva.mmi 2m.fasta > 2m-alignment.paf"

./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_NR99_tax_silva.mmi 1k.fasta > 1k-nr99-alignment.paf"
./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_NR99_tax_silva.mmi 10k.fasta > 10k-nr99-alignment.paf"
./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_NR99_tax_silva.mmi 50k.fasta > 50k-nr99-alignment.paf"
./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_NR99_tax_silva.mmi 100k.fasta > 100k-nr99-alignment.paf"
./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_NR99_tax_silva.mmi 500k.fasta > 500k-nr99-alignment.paf"
#./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_NR99_tax_silva.mmi 1m.fasta > 1m-nr99-alignment.paf"
#./benchmark.sh "minimap2/minimap2" "-x map-ont SILVA_138_SSURef_NR99_tax_silva.mmi 2m.fasta > 2m-nr99-alignment.paf"