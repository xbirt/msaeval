#!/bin/bash

set -e

THREADS=$(nproc)

# Create the kraken2 SILVA database if needed
if [ ! -d "silva_db" ]; then
    export KRAKEN2_DB_NAME=silva_db
    export KRAKEN2_THREAD_CT=32
    cp kraken2/bin/16S_silva_installation.sh . && \
    sed -i 's/SILVA_VERSION=".*"/SILVA_VERSION="138"/' 16S_silva_installation.sh && \
    sed -i 's/SSURef_NR99_tax_silva/SSURef_tax_silva/' 16S_silva_installation.sh && \
    sed -i 's/--threads \$KRAKEN2_THREAD_CT/& --kmer-len 31 --minimizer-len 21 --minimizer-spaces 3/' 16S_silva_installation.sh
    ./16S_silva_installation.sh
fi

./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 1k_output.txt --report 1k_report.txt 1k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 10k_output.txt --report 10k_report.txt 10k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 50k_output.txt --report 50k_report.txt 50k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 100k_output.txt --report 100k_report.txt 100k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 500k_output.txt --report 500k_report.txt 500k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 1m_output.txt --report 1m_report.txt 1m.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 2m_output.txt --report 2m_report.txt 2m.fasta"