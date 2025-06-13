#!/bin/bash

THREADS=$(nproc)

./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 0.1k.fasta > 0.1k-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 0.5k.fasta > 0.5k-alignment.fasta"
./benchmark.sh "clustalo" "-i 0.1k.fasta -o 0.1k-alignment-16.fasta --auto -v --threads=$THREADS --seqtype=RNA --outfmt=fasta"
./benchmark.sh "clustalo" "-i 0.5k.fasta -o 0.5k-alignment-16.fasta --auto -v --threads=$THREADS --seqtype=RNA --outfmt=fasta"
./benchmark.sh "t_coffee" "-thread $THREADS -seq 0.1k.fasta -output=fasta_aln -outfile=0.1k-alignment.fasta"
./benchmark.sh "t_coffee" "-thread $THREADS -seq 0.5k.fasta -output=fasta_aln -outfile=0.5k-alignment.fasta"
./benchmark.sh "/usr/bin/muscle" "-super5 0.1k.fasta -output 0.1k-alignment.fasta"
./benchmark.sh "/usr/bin/muscle" "-super5 0.5k.fasta -output 0.5k-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 0.1k.fasta > 0.1k-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --auto 0.5k.fasta > 0.5k-alignment.fasta"
./benchmark.sh "run_pasta.py" "--num-cpus=$THREADS -d rna -i 0.1k.fasta"
./benchmark.sh "run_pasta.py" "--num-cpus=$THREADS -d rna -i 0.5k.fasta"
./benchmark.sh "kalign" "-n $THREADS --type rna -i 0.1k.fasta -o 0.1k-alignment.fasta"
./benchmark.sh "kalign" "-n $THREADS --type rna -i 0.5k.fasta -o 0.5k-alignment.fasta"
