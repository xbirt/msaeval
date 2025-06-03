#!/bin/bash

set -e

seqkit seq --rna2dna 1k.fasta -o 1k_dna.fasta
mv 1k_dna.fasta 1k.fasta
seqkit seq --rna2dna 10k.fasta -o 10k_dna.fasta
mv 10k_dna.fasta 10k.fasta
seqkit seq --rna2dna 50k.fasta -o 50k_dna.fasta
mv 50k_dna.fasta 50k.fasta
seqkit seq --rna2dna 100k.fasta -o 100k_dna.fasta
mv 100k_dna.fasta 100k.fasta
seqkit seq --rna2dna 500k.fasta -o 500k_dna.fasta
mv 500k_dna.fasta 500k.fasta
seqkit seq --rna2dna 1m.fasta -o 1m_dna.fasta
mv 1m_dna.fasta 1m.fasta
seqkit seq --rna2dna 2m.fasta -o 2m_dna.fasta
mv 2m_dna.fasta 2m.fasta

THREADS=$(nproc)

./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 1k_nr_output.txt --report 1k_nr_report.txt 1k.fasta --classified-out 1k_nr_classified.txt --unclassified-out 1k_nr_unclassified.txt -A silva_nr_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 10k_nr_output.txt --report 10k_nr_report.txt 10k.fasta --classified-out 10k_nr_classified.txt --unclassified-out 10k_nr_unclassified.txt -A silva_nr_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 50k_nr_output.txt --report 50k_nr_report.txt 50k.fasta --classified-out 50k_nr_classified.txt --unclassified-out 50k_nr_unclassified.txt -A silva_nr_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 100k_nr_output.txt --report 100k_nr_report.txt 100k.fasta --classified-out 100k_nr_classified.txt --unclassified-out 100k_nr_unclassified.txt -A silva_nr_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 500k_nr_output.txt --report 500k_nr_report.txt 500k.fasta --classified-out 500k_nr_classified.txt --unclassified-out 500k_nr_unclassified.txt -A silva_nr_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 1m_nr_output.txt --report 1m_nr_report.txt 1m.fasta --classified-out 1m_nr_classified.txt --unclassified-out 1m_nr_unclassified.txt -A silva_nr_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 2m_nr_output.txt --report 2m_nr_report.txt 2m.fasta --classified-out 2m_nr_classified.txt --unclassified-out 2m_nr_unclassified.txt -A silva_nr_db/accession_map.k2map"

./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 1k_output.txt --report 1k_report.txt 1k.fasta --classified-out 1k_classified.txt --unclassified-out 1k_unclassified.txt -A silva_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 10k_output.txt --report 10k_report.txt 10k.fasta --classified-out 10k_classified.txt --unclassified-out 10k_unclassified.txt -A silva_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 50k_output.txt --report 50k_report.txt 50k.fasta --classified-out 50k_classified.txt --unclassified-out 50k_unclassified.txt -A silva_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 100k_output.txt --report 100k_report.txt 100k.fasta --classified-out 100k_classified.txt --unclassified-out 100k_unclassified.txt -A silva_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 500k_output.txt --report 500k_report.txt 500k.fasta --classified-out 500k_classified.txt --unclassified-out 500k_unclassified.txt -A silva_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 1m_output.txt --report 1m_report.txt 1m.fasta --classified-out 1m_classified.txt --unclassified-out 1m_unclassified.txt -A silva_db/accession_map.k2map"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 2m_output.txt --report 2m_report.txt 2m.fasta --classified-out 2m_classified.txt --unclassified-out 2m_unclassified.txt -A silva_db/accession_map.k2map"

# Without accession map
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 1k_nr_output.txt 1k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 10k_nr_output.txt 10k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 50k_nr_output.txt 50k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 100k_nr_output.txt 100k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 500k_nr_output.txt 500k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 1m_nr_output.txt 1m.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_nr_db --threads $(nproc) --output 2m_nr_output.txt 2m.fasta"

./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 1k_output.txt 1k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 10k_output.txt 10k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 50k_output.txt 50k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 100k_output.txt 100k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 500k_output.txt 500k.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 1m_output.txt 1m.fasta"
./benchmark.sh "kraken2/bin/kraken2" "--db silva_db --threads $(nproc) --output 2m_output.txt 2m.fasta"