#!/bin/bash

THREADS=$(nproc)

# Download the SILVA 138 SSU Ref databases
if [ ! -f "SILVA_138_SSURef_NR99_05_01_20_opt.arb" ]; then
    if [ ! -f "SILVA_138_SSURef_NR99_05_01_20_opt.arb.gz" ]; then
        wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/ARB_files/SILVA_138_SSURef_NR99_05_01_20_opt.arb.gz
    fi
    gunzip SILVA_138_SSURef_NR99_05_01_20_opt.arb.gz
fi

if [ ! -f "SILVA_138_SSURef_05_01_20_opt.arb" ]; then
    if [ ! -f "SILVA_138_SSURef_05_01_20_opt.arb.gz" ]; then
        wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/ARB_files/SILVA_138_SSURef_05_01_20_opt.arb.gz
    fi
    gunzip SILVA_138_SSURef_05_01_20_opt.arb.gz
fi

# Pre-build the index for the arb files
if [ ! -f "SILVA_138_SSURef_NR99_05_01_20_opt.sidx" ]; then
    sina/sina --threads 32 --in 0.1k.fasta --db SILVA_138_SSURef_NR99_05_01_20_opt.arb >/dev/null
fi

if [ ! -f "SILVA_138_SSURef_05_01_20_opt.sidx" ]; then
    sina/sina --threads 32 --in 0.1k.fasta --db SILVA_138_SSURef_05_01_20_opt.arb >/dev/null
fi

./benchmark.sh "sina/sina" "--threads $(nproc) --out 1k-align-nr99.fasta --out 1k-nr99.arb --out 1k-nr99.csv --in 1k.fasta --db SILVA_138_SSURef_NR99_05_01_20_opt.arb --log-file 1k.log --meta-fmt csv"
gzip 1k-align-nr99.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 10k-align-nr99.fasta --out 10k-nr99.arb --out 10k-nr99.csv --in 10k.fasta --db SILVA_138_SSURef_NR99_05_01_20_opt.arb --log-file 10k-nr99.log --meta-fmt csv"
gzip 10k-align-nr99.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 50k-align-nr99.fasta --out 50k-nr99.arb --out 50k-nr99.csv --in 50k.fasta --db SILVA_138_SSURef_NR99_05_01_20_opt.arb --log-file 50k-nr99.log --meta-fmt csv"
gzip 50k-align-nr99.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 100k-align-nr99.fasta --out 100k-nr99.arb --out 100k-nr99.csv --in 100k.fasta --db SILVA_138_SSURef_NR99_05_01_20_opt.arb --log-file 100k-nr99.log --meta-fmt csv"
gzip 100k-align-nr99.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 500k-align-nr99.fasta --out 500k-nr99.arb --out 500k-nr99.csv --in 500k.fasta --db SILVA_138_SSURef_NR99_05_01_20_opt.arb --log-file 500k-nr99.log --meta-fmt csv"
gzip 500k-align-nr99.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 1m-align-nr99.fasta --out 1m-nr99.arb --out 1m-nr99.csv --in 1m.fasta --db SILVA_138_SSURef_NR99_05_01_20_opt.arb --log-file 1m-nr99.log --meta-fmt csv"
gzip 1m-align-nr99.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 2m-align-nr99.fasta --out 2m-nr99.arb --out 2m-nr99.csv --in 2m.fasta --db SILVA_138_SSURef_NR99_05_01_20_opt.arb --log-file 2m-nr99.log --meta-fmt csv"
gzip 2m-align-nr99.fasta

./benchmark.sh "sina/sina" "--threads $(nproc) --out 1k-align.fasta --out 1k.arb --out 1k.csv --in 1k.fasta --db SILVA_138_SSURef_05_01_20_opt.arb --log-file 1k.log --meta-fmt csv"
gzip 1k-align.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 10k-align.fasta --out 10k.arb --out 10k.csv --in 10k.fasta --db SILVA_138_SSURef_05_01_20_opt.arb --log-file 10k.log --meta-fmt csv"
gzip 10k-align.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 50k-align.fasta --out 50k.arb --out 50k.csv --in 50k.fasta --db SILVA_138_SSURef_05_01_20_opt.arb --log-file 50k.log --meta-fmt csv"
gzip 50k-align.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 100k-align.fasta --out 100k.arb --out 100k.csv --in 100k.fasta --db SILVA_138_SSURef_05_01_20_opt.arb --log-file 100k.log --meta-fmt csv"
gzip 100k-align.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 500k-align.fasta --out 500k.arb --out 500k.csv --in 500k.fasta --db SILVA_138_SSURef_05_01_20_opt.arb --log-file 500k.log --meta-fmt csv"
gzip 500k-align.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 1m-align.fasta --out 1m.arb --out 1m.csv --in 1m.fasta --db SILVA_138_SSURef_05_01_20_opt.arb --log-file 1m.log --meta-fmt csv"
gzip 1m-align.fasta
./benchmark.sh "sina/sina" "--threads $(nproc) --out 2m-align.fasta --out 2m.arb --out 2m.csv --in 2m.fasta --db SILVA_138_SSURef_05_01_20_opt.arb --log-file 2m.log --meta-fmt csv"
gzip 2m-align.fasta
