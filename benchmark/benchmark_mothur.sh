#!/bin/bash

# Mothur args
# processors - Number of CPU threads to use
# search - Template search method - kmer (default) or suffix
# ksize - Size of kmers for kmer search - 8 (default), range 5-12
# align - Alignment algorithm - needleman (default) or gotoh
# flip - Try reverse complement if alignment is poor - t (true, default) or f (false)
# threshold - Cutoff for flipping sequence - 0.75 (default)


THREADS=$(nproc)
mkdir -p seed_db
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
# The alignment files need to be compressed to save space.
gzip 1k.align
mv 1k.align.gz seed_db/
mv mothur*.logfile seed_db/
mv *.align_report seed_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=10k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
gzip 10k.align
mv 10k.align.gz seed_db/
mv mothur*.logfile seed_db/
mv *.align_report seed_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=50k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
gzip 50k.align
mv 50k.align.gz seed_db/
mv mothur*.logfile seed_db/
mv *.align_report seed_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=100k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
gzip 100k.align
mv 100k.align.gz seed_db/
mv mothur*.logfile seed_db/
mv *.align_report seed_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=500k.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
gzip 500k.align
mv 500k.align.gz seed_db/
mv mothur*.logfile seed_db/
mv *.align_report seed_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1m.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
gzip 1m.align
mv 1m.align.gz seed_db/
mv mothur*.logfile seed_db/
mv *.align_report seed_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=2m.fasta, reference=silva.seed_v138_2.align, processors=$THREADS)\""
gzip 2m.align
mv 2m.align.gz seed_db/
mv mothur*.logfile seed_db/
mv *.align_report seed_db/

mkdir -p nr_db
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1k.fasta, reference=silva.nr_v138_2.align, processors=$THREADS)\""
gzip 1k.align
mv 1k.align.gz nr_db/
mv mothur*.logfile nr_db/
mv *.align_report nr_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=10k.fasta, reference=silva.nr_v138_2.align, processors=$THREADS)\""
gzip 10k.align
mv 10k.align.gz nr_db/
mv mothur*.logfile nr_db/
mv *.align_report nr_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=50k.fasta, reference=silva.nr_v138_2.align, processors=$THREADS)\""
gzip 50k.align
mv 50k.align.gz nr_db/
mv mothur*.logfile nr_db/
mv *.align_report nr_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=100k.fasta, reference=silva.nr_v138_2.align, processors=$THREADS)\""
gzip 100k.align
mv 100k.align.gz nr_db/
mv mothur*.logfile nr_db/
mv *.align_report nr_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=500k.fasta, reference=silva.nr_v138_2.align, processors=$THREADS)\""
gzip 500k.align
mv 500k.align.gz nr_db/
mv mothur*.logfile nr_db/
mv *.align_report nr_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1m.fasta, reference=silva.nr_v138_2.align, processors=$THREADS)\""
gzip 1m.align
mv 1m.align.gz nr_db/
mv mothur*.logfile nr_db/
mv *.align_report nr_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=2m.fasta, reference=silva.nr_v138_2.align, processors=$THREADS)\""
gzip 2m.align
mv 2m.align.gz nr_db/
mv mothur*.logfile nr_db/
mv *.align_report nr_db/

mkdir -p nr99_db
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1k.fasta, reference=SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 1k.align nr99_db/
mv mothur*.logfile nr99_db/
mv *.align_report nr99_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=10k.fasta, reference=SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 10k.align nr99_db/
mv mothur*.logfile nr99_db/
mv *.align_report nr99_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=50k.fasta, reference=SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 50k.align nr99_db/
mv mothur*.logfile nr99_db/
mv *.align_report nr99_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=100k.fasta, reference=SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 100k.align nr99_db/
mv mothur*.logfile nr99_db/
mv *.align_report nr99_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=500k.fasta, reference=SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 500k.align nr99_db/
mv mothur*.logfile nr99_db/
mv *.align_report nr99_db/
#./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1m.fasta, reference=SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
#mv 1m.align nr99_db/
#mv mothur*.logfile nr99_db/
#mv *.align_report nr99_db/
#./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=2m.fasta, reference=SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
#mv 2m.align nr99_db/
#mv mothur*.logfile nr99_db/
#mv *.align_report nr99_db/

mkdir -p full_db
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1k.fasta, reference=SILVA_138_SSURef_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 1k.align full_db/
mv mothur*.logfile full_db/
mv *.align_report full_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=10k.fasta, reference=SILVA_138_SSURef_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 10k.align full_db/
mv mothur*.logfile full_db/
mv *.align_report full_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=50k.fasta, reference=SILVA_138_SSURef_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 50k.align full_db/
mv mothur*.logfile full_db/
mv *.align_report full_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=100k.fasta, reference=SILVA_138_SSURef_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 100k.align full_db/
mv mothur*.logfile full_db/
mv *.align_report full_db/
./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=500k.fasta, reference=SILVA_138_SSURef_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
mv 500k.align full_db/
mv mothur*.logfile full_db/
mv *.align_report full_db/
#./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=1m.fasta, reference=SILVA_138_SSURef_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
#mv 1m.align full_db/
#mv mothur*.logfile full_db/
#mv *.align_report full_db/
#./benchmark.sh "mothur/mothur" "\"#align.seqs(fasta=2m.fasta, reference=SILVA_138_SSURef_tax_silva_full_align_trunc.fasta, processors=$THREADS)\""
#mv 2m.align full_db/
#mv mothur*.logfile full_db/
#mv *.align_report full_db/
