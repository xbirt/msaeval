#!/bin/bash

# Script pentru generarea de citiri cu grinder și post-procesarea cu CuReSim-LoRM
# Cerințe:
#   - SeqKit instalat și accesibil în PATH
#   - grinder instalat și accesibil în PATH
#   - CuReSim-LoRM în format jar localizat în CuReSim-LoRM/CuReSim-LoRM.jar
#   - parallel instalat

# Aflăm numărul de nuclee ale procesorului
MAX_CORES=$(nproc)

INPUT_FILE="selected_species_filtered.fasta"
SUBSET_FILE="selected_species_subset.fasta"

RANDOM_SEED=42
PRIMER_FILE="primers_27F_1492R.fasta"

echo "Se extrag subseturi din fișierele disponibile..."
# Extragem diferite procente din fiecare fișier
[ ! -f "selected_species_1000_subset.fasta" ] && seqkit sample -p 0.001 -s $RANDOM_SEED "selected_species_1000_filtered.fasta" > "selected_species_1000_subset.fasta"
[ ! -f "selected_species_10000_subset.fasta" ] && seqkit sample -p 0.01 -s $RANDOM_SEED "selected_species_10000_filtered.fasta" > "selected_species_10000_subset.fasta"
[ ! -f "selected_species_50000_subset.fasta" ] && seqkit sample -p 0.05 -s $RANDOM_SEED "selected_species_50000_filtered.fasta" > "selected_species_50000_subset.fasta"
[ ! -f "selected_species_100000_subset.fasta" ] && seqkit sample -p 0.1 -s $RANDOM_SEED "selected_species_100000_filtered.fasta" > "selected_species_100000_subset.fasta"
[ ! -f "selected_species_500000_subset.fasta" ] && seqkit sample -p 0.1 -s $RANDOM_SEED "selected_species_500000_filtered.fasta" > "selected_species_500000_subset.fasta"
[ ! -f "selected_species_1000000_subset.fasta" ] && seqkit sample -p 0.1 -s $RANDOM_SEED "selected_species_1000000_filtered.fasta" > "selected_species_1000000_subset.fasta"

echo "Se crează fișierul cu primers..."
# Folosim primerii V1-V9 folosiți și în kit-ul comercial Nanopore 16S Barcoding Kit 24 V14 (27F and 1492R)
[ ! -f "$PRIMER_FILE" ] && cat > "$PRIMER_FILE" << 'EOF'
>forward_primer;27F
AGAGTTTGATCMTGGCTCAG
>reverse_primer;1492R
TACGGYTACCTTGTTACGACTT
EOF

echo "Se generează citirile cu grinder..."
# Executăm grinder în paralel cu diferite cantități de citiri generate
parallel --jobs $MAX_CORES --link "grinder -reference_file selected_species_{4}_subset.fasta \
  -forward_reverse $PRIMER_FILE \
  -unidirectional 1 \
  -total_reads {1} \
  -abundance_model logarithmic 1.5 \
  -mutation_dist uniform 0.3 \
  -mutation_ratio 90 10 \
  -read_dist 1500 normal 30 \
  -chimera_perc 10 \
  -random_seed {2} \
  -fastq_output 1 \
  -qual_levels 30 10 \
  -base_name {3} > {3}_grinder.log 2>&1" ::: 1000 10000 50000 100000 200000 500000 1000000 2000000 :::+ $(seq $INITIAL_SEED $((INITIAL_SEED+7))) ::: 1k 10k 50k 100k 200k 500k 1m 2m ::: 1000 10000 50000 100000 100000 500000 1000000 1000000

echo "Fișierele cu citiri grinder au fost finalizate."

echo "Se filtrează citirile mai scurte de 1000 pb..."
# Filtrăm citirile scurte (sub 1000 nucleotide) și ștergem sursa
parallel --jobs $MAX_CORES "seqkit seq -m 1000 {} -o {.}-filtered.fastq && rm {}" ::: *-reads.fastq

echo "Se generează fișierul profil pentru CuReSim-LoRM..."
# Generăm profil cu numărul maxim de înregistrări
max_num_reads=2000000
echo "read;ins;del;sub" > nanopore_v14.txt
for i in $(seq 1 "$max_num_reads"); do
    echo "$i;0.004;0.003;0.001" >> nanopore_v14.txt
done
echo "Fișierul profil pentru CuReSim-LoRM a fost finalizat."

# CuReSim-LoRM nu poate fi executat în paralel deoarece folosește fișiere temporare cu nume predefinit.
# Procesare secvențială CuReSim-LoRM
for moniker in 1k 10k 50k 100k 200k 500k 1m 2m; do
    input_file="${moniker}-reads-filtered.fastq"
    output_file="curesim-${moniker}.fastq"
    
    # Dacă fișierul de ieșire există, sărim acest pas
    if [ -f "$output_file" ] && [ -s "$output_file" ]; then
        echo "Sărim peste ${moniker} - fișierul de ieșire există deja."
        continue
    fi
    
    echo "Se execută CuReSim-LoRM pe ${moniker}"
    
    # CuReSim-LoRM are probleme cu fișierele mari, deci împărțim fișierul în bucăți de 100000 de citiri pe care știm că le poate procesa
    echo "Se împarte fișierul de intrare în părți de maxim 100000 de citiri..."
    seqkit split2 -s 100000 "$input_file"
    
    # Procesăm fiecare bucată cu CuReSim-LoRM
    for part_file in "${input_file}.split/"*.fastq; do
        part_num=$(echo "$part_file" | grep -oP 'part_\K\d+')
        output_part="curesim-${moniker}-mix.part_${part_num}.fastq"
        
        echo "Se procesează partea ${part_num}..."
        java -jar CuReSim-LoRM/CuReSim-LoRM.jar \
            -Xms16384m -Xmx24576m \
            -f "$part_file" \
            -p nanopore_v14.txt \
            -o "$output_part" \
            --keep_headers > "${moniker}_curesim_part_${part_num}.log" 2>&1
    done
    
    # Combinăm toate părțile procesate
    echo "Se combină părțile procesate..."
    cat curesim-${moniker}-mix.part_*.fastq > "curesim-${moniker}-mix.fastq"
    
    # Ștergem fișierele partiționate
    rm -f curesim-${moniker}-mix.part_*.fastq
    rm -rf "${input_file}.split"
    
    # CuReSim-LoRM folosește T (ADN) în loc de U deși noi avem ARN în sursă, așadar trebuie să convertim T în U și să formatăm ieșirea
    seqkit seq -w 0 "curesim-${moniker}-mix.fastq" | sed '/^@/!{/^+/!s/T/U/g}' | seqkit seq -w 0 -o "$output_file"
    rm -f "curesim-${moniker}-mix.fastq"
done

# Ștergem fișier profil
rm nanopore_v14.txt
echo "Fișierul profil pentru CuReSim-LoRM a fost șters."

echo "Se filtrează citirile scurte din fișierele curesim și se face conversie FASTQ la FASTA..."
parallel --jobs $MAX_CORES "seqkit seq -m 1000 --max-len 2000 {} -o {.}-filtered.fastq" ::: curesim-*.fastq
parallel --jobs $MAX_CORES "seqkit fq2fa {} -o {.}.fasta" ::: curesim-*-filtered.fastq

echo "Se redenumesc și se compresează fișierele FASTA..."
for file in curesim-*-filtered.fasta; do
    # Extragem X dintre "curesim-" și "-filtered.fasta"
    x=$(echo "$file" | sed -E 's/^curesim-(.*)-filtered\.fasta$/\1/')
    # Redenumim fișierul în X.fasta
    mv "$file" "$x.fasta"
done
find . -maxdepth 1 -regextype posix-extended -regex './[0-9]+[km]\.fasta' -print | tar -cvf reads.tar -T -
bzip2 reads.tar