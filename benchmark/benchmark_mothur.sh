#!/bin/bash

THREADS=$(nproc)

./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=0.1k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=0.5k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=10k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=50k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=100k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=500k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1m.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=2m.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""

# Args
# processors - Number of CPU threads to use
# search - Template search method - kmer (default) or suffix
# ksize - Size of kmers for kmer search - 8 (default), range 5-12
# align - Alignment algorithm - needleman (default) or gotoh
# flip - Try reverse complement if alignment is poor - t (true, default) or f (false)
# threshold - Cutoff for flipping sequence - 0.75 (default)
