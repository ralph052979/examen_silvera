#!/bin/bash
# =============================================================================
#  parteA1_estado.sh - estado de las corridas de la Tarea A.1 (no gasta horas-nucleo).
#  Uso:  bash slurm/parteA1_estado.sh
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
REG=logs/jobs_A1.txt
[ -s "$REG" ] || { echo "Todavia no hay corridas de la A.1 (no existe $REG)."; exit 0; }

printf "%-3s %-3s %-3s %-9s %-11s %-19s %8s %-8s %-11s  %s\n" \
       "#" r p Job Slurm Nodo Minutos Afinidad "Carga ajena" Resultado
I=0; NOK=0; NMAL=0; NCURSO=0
while read -r P R ID F H; do
    [ -n "$P" ] || continue
    I=$((I + 1))
    LINEA=$(squeue -h -j "$ID" -o "%T %M" 2>/dev/null | head -1)
    if [ -n "$LINEA" ]; then
        ST=$(echo "$LINEA" | awk '{print $1}'); MIN=$(echo "$LINEA" | awk '{print $2}')
    else
        ST=$(sacct -n -X -j "$ID" -o State%20 2>/dev/null | head -1 | tr -d ' '); MIN="-"
    fi
    [ -n "$ST" ] || ST="?"
    D="corridas/ex_A1_p${P}_r${R}_${ID}"
    NODO="-"; AF="-"; AJ="-"; RES="-"
    if [ -d "$D" ]; then
        NODO=$(grep '^nodo=' "$D/entorno.txt" 2>/dev/null | cut -d= -f2 | cut -d. -f1)
        if [ -s "$D/tiempo.txt" ]; then
            W=$(grep '^wall_s=' "$D/tiempo.txt" | cut -d= -f2)
            MIN=$(awk -v w="${W:-0}" 'BEGIN{printf "%.1f", w/60}')
            AF=$(grep '^afinidad_ok=' "$D/tiempo.txt" | cut -d= -f2)
            [ "$AF" = "1" ] && AF="ok" || AF="MAL"
            AJ=$(grep '^carga_ajena_max=' "$D/tiempo.txt" | cut -d= -f2)
            RES=$(grep '^estado=' "$D/tiempo.txt" | tail -1 | cut -d= -f2)
        fi
    fi
    case "$ST" in
        PENDING|RUNNING|CONFIGURING|COMPLETING) NCURSO=$((NCURSO + 1)) ;;
        *) if [ "$RES" = "OK" ]; then NOK=$((NOK + 1)); else NMAL=$((NMAL + 1)); [ "$RES" = "-" ] && RES="$ST"; fi ;;
    esac
    printf "%-3s %-3s %-3s %-9s %-11s %-19s %8s %-8s %-11s  %s\n" \
           "$I" "r$R" "$P" "$ID" "$ST" "${NODO:--}" "$MIN" "$AF" "${AJ:--}" "$RES"
done < "$REG"
echo
echo "Terminadas bien: $NOK   con problema: $NMAL   en cola o corriendo: $NCURSO"
if [ "$NCURSO" -eq 0 ]; then
    echo "Todas terminaron. Siguiente paso:  python3 analisis/parteA1_resumen.py"
else
    echo "La cadena corre una a la vez: las que dicen PENDING esperan a la anterior (es normal)."
fi
[ "$NMAL" -eq 0 ] || echo "Hay corridas con problema: mira logs/ex_A1_p<p>_r<r>_<Job>.out y mandame una captura."
