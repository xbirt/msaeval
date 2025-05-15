#!/bin/bash

THREADS=$(nproc)

./benchmark.sh "run_pasta.py" "--num-cpus=$THREADS -d rna -i 1k.fasta"
./benchmark.sh "run_pasta.py" "--num-cpus=$THREADS -d rna -i 10k.fasta"
#./benchmark.sh "run_pasta.py" "--num-cpus=$THREADS -d rna -i 50k.fasta"
#./benchmark.sh "run_pasta.py" "--num-cpus=$THREADS -d rna -i 100k.fasta"
#./benchmark.sh "run_pasta.py" "--num-cpus=$THREADS -d rna -i 500k.fasta"
#./benchmark.sh "run_pasta.py" "--num-cpus=$THREADS -d rna -i 1m.fasta"
#./benchmark.sh "run_pasta.py" "--num-cpus=$THREADS -d rna -i 2m.fasta"
