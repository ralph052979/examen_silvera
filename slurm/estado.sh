#!/bin/bash
# =============================================================================
#  estado.sh - estado de los jobs de la Tarea 0 (no gasta horas-nucleo).
#  Uso:  bash slurm/estado.sh
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
REG=logs/jobs_t0.txt
[ -s "$REG" ] || { echo "Todavia no hay jobs enviados (no existe $REG)."; exit 0; }

printf "%-4s %-9s %-11s %-20s %5s %7s  %s\n" Cand Job Slurm Nodo Iter Minutos Resultado
NOK=0; NMAL=0; NCURSO=0; NTIME=0
while read -r K ID F H; do
    [ -n "$K" ] || continue
    LINEA=$(squeue -h -j "$ID" -o "%T %M %N" 2>/dev/null | head -1)
    if [ -n "$LINEA" ]; then
        ST=$(echo "$LINEA" | awk '{print $1}')
        MIN=$(echo "$LINEA" | awk '{print $2}')
    else
        ST=$(sacct -n -X -j "$ID" -o State%20 2>/dev/null | head -1 | tr -d ' ')
        MIN="-"
    fi
    [ -n "$ST" ] || ST="?"
    D="corridas/ex_T0b_${K}_${ID}"
    NODO="-"; IT="-"; RES="-"
    if [ -d "$D" ]; then
        NODO=$(grep '^nodo=' "$D/entorno.txt" 2>/dev/null | cut -d= -f2)
        IT=$(grep -cE '^(DAV|RMM):' "$D/OSZICAR" 2>/dev/null)
        IT=${IT:-0}
        if [ -s "$D/tiempo.txt" ]; then
            W=$(grep '^wall_s=' "$D/tiempo.txt" | cut -d= -f2)
            MIN=$(awk -v w="${W:-0}" 'BEGIN{printf "%.1f", w/60}')
        fi
        RES=$(grep '^estado=' "$D/tiempo.txt" 2>/dev/null | tail -1 | cut -d= -f2)
        RES=${RES:--}
    fi
    case "$ST" in
        PENDING|RUNNING|CONFIGURING|COMPLETING) NCURSO=$((NCURSO + 1)) ;;
        TIMEOUT) NTIME=$((NTIME + 1)); RES="se paso del limite de tiempo (t1 > 80 min)" ;;
        *) if [ "$RES" = "OK" ]; then NOK=$((NOK + 1)); else NMAL=$((NMAL + 1)); fi ;;
    esac
    printf "%-4s %-9s %-11s %-20s %5s %7s  %s\n" "$K" "$ID" "$ST" "${NODO:--}" "$IT" "$MIN" "$RES"
done < "$REG"
echo
echo "Terminados bien: $NOK   con problema: $NMAL   en cola o corriendo: $NCURSO"
[ "$NTIME" -eq 0 ] || echo "Se pasaron del limite de tiempo: $NTIME (no es un error: esa malla es demasiado lenta para la ventana)"
if [ "$NCURSO" -eq 0 ]; then
    echo "Todos terminaron. Siguiente paso:  python3 analisis/resumen_t0.py"
fi
[ "$NMAL" -eq 0 ] || echo "Hay jobs con problema: mira logs/ex_T0b_<Cand>_<Job>.out y mandame una captura."
