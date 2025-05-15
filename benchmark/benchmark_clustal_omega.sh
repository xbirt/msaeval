#!/bin/bash

THREADS=$(nproc)

./benchmark.sh "clustalo" "-i 1k.fasta -o 1k-alignment.fasta --auto -v --threads=$THREADS --seqtype=RNA --outfmt=fasta"
./benchmark.sh "clustalo" "-i 10k.fasta -o 10k-alignment.fasta --auto -v --threads=$THREADS --seqtype=RNA --outfmt=fasta"
./benchmark.sh "clustalo" "-i 50k.fasta -o 50k-alignment.fasta --auto -v --threads=$THREADS --seqtype=RNA --outfmt=fasta"
./benchmark.sh "clustalo" "-i 100k.fasta -o 100k-alignment.fasta --auto -v --threads=$THREADS --seqtype=RNA --outfmt=fasta"
#./benchmark.sh "clustalo" "-i 500k.fasta -o 500k-alignment.fasta --auto -v --threads=$THREADS --seqtype=RNA --outfmt=fasta"
#./benchmark.sh "clustalo" "-i 1m.fasta -o 1m-alignment.fasta --auto -v --threads=$THREADS --seqtype=RNA --outfmt=fasta"
#./benchmark.sh "clustalo" "-i 2m.fasta -o 2m-alignment.fasta --auto -v --threads=$THREADS --seqtype=RNA --outfmt=fasta"
