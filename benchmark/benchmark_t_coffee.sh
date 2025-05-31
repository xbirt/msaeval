#!/bin/bash

export MAX_N_PID_4_TCOFFEE=$(cat /proc/sys/kernel/pid_max)

THREADS=$(nproc)

./benchmark.sh "t_coffee" "-cache=no -thread $THREADS -seq 1k.fasta -output=fasta_aln -outfile=1k-alignment.fasta"
./benchmark.sh "t_coffee" "-cache=no -thread $THREADS -seq 10k.fasta -output=fasta_aln -outfile=10k-alignment.fasta"
./benchmark.sh "t_coffee" "-cache=no -thread $THREADS -seq 50k.fasta -output=fasta_aln -outfile=50k-alignment.fasta"
./benchmark.sh "t_coffee" "-cache=no -thread $THREADS -seq 100k.fasta -output=fasta_aln -outfile=100k-alignment.fasta"
#./benchmark.sh "t_coffee" "-cache=no -thread $THREADS -seq 500k.fasta -output=fasta_aln -outfile=500k-alignment.fasta"
#./benchmark.sh "t_coffee" "-cache=no -thread $THREADS -seq 1m.fasta -output=fasta_aln -outfile=1m-alignment.fasta"
#./benchmark.sh "t_coffee" "-cache=no -thread $THREADS -seq 2m.fasta -output=fasta_aln -outfile=2m-alignment.fasta"
