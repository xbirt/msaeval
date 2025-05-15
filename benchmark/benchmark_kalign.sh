#!/bin/bash

THREADS=$(nproc)

./benchmark.sh "kalign" "-n $THREADS --type rna -i 1k.fasta -o 1k-alignment.fasta"
./benchmark.sh "kalign" "-n $THREADS --type rna -i 10k.fasta -o 10k-alignment.fasta"
./benchmark.sh "kalign" "-n $THREADS --type rna -i 50k.fasta -o 50k-alignment.fasta"
#./benchmark.sh "kalign" "-n $THREADS --type rna -i 100k.fasta -o 100k-alignment.fasta"
#./benchmark.sh "kalign" "-n $THREADS --type rna -i 500k.fasta -o 500k-alignment.fasta"
#./benchmark.sh "kalign" "-n $THREADS --type rna -i 1m.fasta -o 1m-alignment.fasta"
#./benchmark.sh "kalign" "-n $THREADS --type rna -i 2m.fasta -o 2m-alignment.fasta"
