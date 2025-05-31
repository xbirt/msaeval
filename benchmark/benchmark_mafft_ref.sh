#!/bin/bash

THREADS=$(nproc)

# Replace dots with dashes in the SILVA databases. Mafft is not able to handle dots in the sequences.
#sed 's/\./-/g' silva.seed_v138_2.align > silva.seed_v138_2_mafft.align
#sed 's/\./-/g' silva.nr_v138_2.align > silva.nr_v138_2_mafft.align

#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 1k.fasta silva.seed_v138_2_mafft.align > 1k-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 10k.fasta silva.seed_v138_2_mafft.align > 10k-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 50k.fasta silva.seed_v138_2_mafft.align > 50k-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 100k.fasta silva.seed_v138_2_mafft.align > 100k-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 500k.fasta silva.seed_v138_2_mafft.align > 500k-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 1m.fasta silva.seed_v138_2_mafft.align > 1m-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 2m.fasta silva.seed_v138_2_mafft.align > 2m-ref-alignment.fasta"

./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 1k.fasta silva.nr_v138_2_mafft.align > 1k-ref-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 10k.fasta silva.nr_v138_2_mafft.align > 10k-ref-alignment.fasta"
./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 50k.fasta silva.nr_v138_2_mafft.align > 50k-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 100k.fasta silva.nr_v138_2_mafft.align > 100k-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 500k.fasta silva.nr_v138_2_mafft.align > 500k-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 1m.fasta silva.nr_v138_2_mafft.align > 1m-ref-alignment.fasta"
#./benchmark.sh "mafft" "--thread $THREADS --quiet --addfragments 2m.fasta silva.nr_v138_2_mafft.align > 2m-ref-alignment.fasta"
