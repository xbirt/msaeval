#!/bin/bash

./benchmark.sh "/usr/bin/muscle" "-super5 1k.fasta -output 1k-alignment.fasta"
./benchmark.sh "/usr/bin/muscle" "-super5 10k.fasta -output 10k-alignment.fasta"
./benchmark.sh "/usr/bin/muscle" "-super5 50k.fasta -output 50k-alignment.fasta"
./benchmark.sh "/usr/bin/muscle" "-super5 100k.fasta -output 100k-alignment.fasta"
#./benchmark.sh "/usr/bin/muscle" "-super5 500k.fasta -output 500k-alignment.fasta"
#./benchmark.sh "/usr/bin/muscle" "-super5 1m.fasta -output 1m-alignment.fasta"
#./benchmark.sh "/usr/bin/muscle" "-super5 2m.fasta -output 2m-alignment.fasta"

# Align is slower, but more accurate
#./benchmark.sh "/usr/bin/muscle" "-align 1k.fasta -output 1k-alignment.fasta"
#./benchmark.sh "/usr/bin/muscle" "-align 10k.fasta -output 10k-alignment.fasta"
#./benchmark.sh "/usr/bin/muscle" "-align 50k.fasta -output 50k-alignment.fasta"
#./benchmark.sh "/usr/bin/muscle" "-align 100k.fasta -output 100k-alignment.fasta"
#./benchmark.sh "/usr/bin/muscle" "-align 500k.fasta -output 500k-alignment.fasta"
#./benchmark.sh "/usr/bin/muscle" "-align 1m.fasta -output 1m-alignment.fasta"
#./benchmark.sh "/usr/bin/muscle" "-align 2m.fasta -output 2m-alignment.fasta"
