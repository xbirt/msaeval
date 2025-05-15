#!/bin/bash

# Utilizare: ./benchmark.sh <unealtă_aliniere> [argumente...]
# Dacă comanda necesită redirecționare output: ./benchmark.sh "<unealtă_aliniere>" "[argumente...] > output_file"
# Exemplu: ./benchmark.sh "mafft" "--thread 32 --quiet --auto 1k.fasta > 1k-alignment.fasta"

TOOL="$1"
shift
ARGUMENTE=("$@")  # Folosim array pentru argumente

OUT_PREFIX="benchmark_$(date +%s)"
LOCK_FILE="/tmp/${OUT_PREFIX}_lock.tmp"
JURNAL_TIMP="${OUT_PREFIX}_timp.log"
JURNAL_DISC="${OUT_PREFIX}_disc.log"
JURNAL_DETALII="${OUT_PREFIX}_detalii.log"
IESIRE_NORMALA="${OUT_PREFIX}_iesire.log"
IESIRE_ERORI="${OUT_PREFIX}_erori.log"

# Inițializăm jurnalul de detalii
echo "Comanda: $TOOL ${ARGUMENTE[*]}" > "$JURNAL_DETALII"
echo "Număr de nuclee: $(nproc)" >> "$JURNAL_DETALII"
echo "Memorie RAM (GB): $(free -g | awk '/^Mem:/{print $2}')" >> "$JURNAL_DETALII"
echo "Memorie RAM (MB): $(free -m | awk '/^Mem:/{print $2}')" >> "$JURNAL_DETALII"
echo "Ora start: $(date)" >> "$JURNAL_DETALII"
echo "Comanda: $TOOL ${ARGUMENTE[*]}"
echo "Ora start: $(date)"

# Creare fișier lock
touch "$LOCK_FILE"

# Măsurăm spațiul inițial pe disc
SPATIU_INITIAL_ROOT=$(df -k / | awk 'NR==2 {print $3}')
SPATIU_INITIAL_TMP=$(df -k /tmp | awk 'NR==2 {print $3}')

# Pornim monitorizarea discului în subshell
(
    MAX_ROOT=$SPATIU_INITIAL_ROOT
    MAX_TMP=$SPATIU_INITIAL_TMP
    
    while [ -f "$LOCK_FILE" ]; do
        CUR_ROOT=$(df -k / | awk 'NR==2 {print $3}')
        CUR_TMP=$(df -k /tmp | awk 'NR==2 {print $3}')
        
        ((CUR_ROOT > MAX_ROOT)) && MAX_ROOT=$CUR_ROOT
        ((CUR_TMP > MAX_TMP)) && MAX_TMP=$CUR_TMP
        
        sleep 0.2
    done
    
    # Scriere rezultate
    {
        echo "SPATIU_MAX_ROOT: $MAX_ROOT"
        echo "SPATIU_MAX_TMP: $MAX_TMP" 
        echo "SPATIU_INITIAL_ROOT: $SPATIU_INITIAL_ROOT"
        echo "SPATIU_INITIAL_TMP: $SPATIU_INITIAL_TMP"
    } > "$JURNAL_DISC"
) &
PID_MONITOR=$!

# Verificăm dacă comanda conține redirecționări
if [[ "$ARGUMENTE" =~ [\>\|] ]]; then
    eval "/usr/bin/time -v -o \"$JURNAL_TIMP\" \"$TOOL\" $ARGUMENTE 2> \"$IESIRE_ERORI\""
else
    eval "/usr/bin/time -v -o \"$JURNAL_TIMP\" \"$TOOL\" $ARGUMENTE > \"$IESIRE_NORMALA\" 2> \"$IESIRE_ERORI\""
fi
COD_RETURN=$?

# Măsurăm spațiul final pe disc
SPATIU_FINAL_ROOT=$(df -k / | awk 'NR==2 {print $3}')
SPATIU_FINAL_TMP=$(df -k /tmp | awk 'NR==2 {print $3}')

# Oprim monitorizarea discului
rm -f "$LOCK_FILE"
wait $PID_MONITOR 2>/dev/null

# Procesăm rezultatele
if [[ -f "$JURNAL_DISC" ]]; then
    DELTA_ROOT=$(( $(grep SPATIU_MAX_ROOT "$JURNAL_DISC" | awk '{print $2}') - SPATIU_INITIAL_ROOT ))
    DELTA_TMP=$(( $(grep SPATIU_MAX_TMP "$JURNAL_DISC" | awk '{print $2}') - SPATIU_INITIAL_TMP ))
    DELTA_FINAL_ROOT=$(( $(grep SPATIU_MAX_ROOT "$JURNAL_DISC" | awk '{print $2}') - SPATIU_FINAL_ROOT ))
    DELTA_FINAL_TMP=$(( $(grep SPATIU_MAX_TMP "$JURNAL_DISC" | awk '{print $2}') - SPATIU_FINAL_TMP ))
    
    echo "=== Rezultate benchmark ==="
    grep -E "Elapsed \(wall clock\) time|User time|System time|Maximum resident set size" "$JURNAL_TIMP"
    
    echo -e "\nSpațiu disc consumat:"
    echo "  - /:    $DELTA_FINAL_ROOT KB"
    echo "  - /tmp: $DELTA_FINAL_TMP KB"

    # Adăugăm delta în jurnalul de disc
    echo "DELTA_ROOT: $DELTA_ROOT KB" >> "$JURNAL_DISC"
    echo "DELTA_TMP: $DELTA_TMP KB" >> "$JURNAL_DISC"
    echo "DELTA_FINAL_ROOT: $DELTA_FINAL_ROOT KB" >> "$JURNAL_DISC"
    echo "DELTA_FINAL_TMP: $DELTA_FINAL_TMP KB" >> "$JURNAL_DISC"
else
    echo "Eroare: Monitorizarea discului a eșuat!" >&2
fi

# Finalizăm jurnalul
echo "Ora sfârșit: $(date)" >> "$JURNAL_DETALII"
echo "Cod ieșire: $COD_RETURN" >> "$JURNAL_DETALII"
echo "Ora sfârșit: $(date)"

exit $COD_RETURN
