#!/bin/bash

THREADS=$(nproc)

./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 1k.fasta > 1k-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 10k.fasta > 10k-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 50k.fasta > 50k-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 100k.fasta > 100k-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 500k.fasta > 500k-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 1m.fasta > 1m-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 2m.fasta > 2m-alignment.fasta"
