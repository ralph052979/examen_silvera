#!/bin/bash
# =============================================================================
#  parteA2_estado.sh - estado del diagnostico y del barrido de la Tarea A.2
#                      (no gasta horas-nucleo).
#  Uso:  bash slurm/parteA2_estado.sh
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

estado_job() {   # imprime "ESTADO MINUTOS"
    local L
    L=$(squeue -h -j "$1" -o "%T %M" 2>/dev/null | head -1)
    if [ -n "$L" ]; then echo "$L"; else echo "$(sacct -n -X -j "$1" -o State%20 2>/dev/null | head -1 | tr -d ' ') -"; fi
}

if [ -s logs/jobs_A2d.txt ]; then
    echo "=== Diagnostico de OpenMP ==="
    printf "%-5s %-9s %-11s %-12s %8s %-9s %-7s  %s\n" Var Job Slurm Nodo Minutos Mascara Hilos Resultado
    NC=0
    while read -r V P ID RESTO; do
        [ -n "$V" ] || continue
        if [ "$ID" = "RECHAZADO" ]; then
            printf "%-5s %-9s %-11s %-12s %8s %-9s %-7s  %s\n" "$V$P" - RECHAZADO - - - - "$(echo "$RESTO" | cut -c1-60)"
            continue
        fi
        read -r ST MIN <<< "$(estado_job "$ID")"
        [ -n "$ST" ] || ST="?"
        case "$ST" in PENDING|RUNNING|CONFIGURING|COMPLETING) NC=$((NC + 1)) ;; esac
        D="corridas/ex_A2d_${V}${P}_${ID}"
        NODO="-"; MAS="-"; HIL="-"; RES="-"
        if [ -s "$D/tiempo.txt" ]; then
            NODO=$(grep '^nodo=' "$D/entorno.txt" 2>/dev/null | cut -d= -f2 | cut -d. -f1 | sed 's/granado-//')
            W=$(grep '^wall_s=' "$D/tiempo.txt" | cut -d= -f2)
            MIN=$(awk -v w="${W:-0}" 'BEGIN{printf "%.1f", w/60}')
            MAS="$(grep '^nucleos_fisicos_distintos=' "$D/tiempo.txt" | cut -d= -f2)/$P"
            HM=$(grep '^hilos_muestras=' "$D/tiempo.txt" | cut -d= -f2)
            [ "${HM:-0}" -gt 0 ] && HIL="$(grep '^nucleos_fisicos_hilos=' "$D/tiempo.txt" | cut -d= -f2)/$P"
            RES=$(grep '^estado=' "$D/tiempo.txt" | tail -1 | cut -d= -f2)
        fi
        printf "%-5s %-9s %-11s %-12s %8s %-9s %-7s  %s\n" "$V$P" "$ID" "$ST" "$NODO" "$MIN" "$MAS" "$HIL" "$RES"
    done < logs/jobs_A2d.txt
    echo
    if [ "$NC" -eq 0 ]; then
        echo "Diagnostico terminado. Siguiente paso:  python3 analisis/parteA2_diagnostico.py"
    else
        echo "Quedan $NC corridas del diagnostico en cola o corriendo."
    fi
    echo
fi

if [ -s logs/jobs_A2.txt ]; then
    echo "=== Barrido OpenMP ==="
    printf "%-3s %-3s %-3s %-9s %-11s %-12s %8s %-7s %-11s  %s\n" "#" r p Job Slurm Nodo Minutos Hilos "Carga ajena" Resultado
    I=0; NOK=0; NMAL=0; NCURSO=0
    while read -r P R ID F H; do
        [ -n "$P" ] || continue
        I=$((I + 1))
        read -r ST MIN <<< "$(estado_job "$ID")"
        [ -n "$ST" ] || ST="?"
        D="corridas/ex_A2_p${P}_r${R}_${ID}"
        NODO="-"; HIL="-"; AJ="-"; RES="-"
        if [ -s "$D/tiempo.txt" ]; then
            NODO=$(grep '^nodo=' "$D/entorno.txt" 2>/dev/null | cut -d= -f2 | cut -d. -f1 | sed 's/granado-//')
            W=$(grep '^wall_s=' "$D/tiempo.txt" | cut -d= -f2)
            MIN=$(awk -v w="${W:-0}" 'BEGIN{printf "%.1f", w/60}')
            HM=$(grep '^hilos_muestras=' "$D/tiempo.txt" | cut -d= -f2)
            [ "${HM:-0}" -gt 0 ] && HIL="$(grep '^nucleos_fisicos_hilos=' "$D/tiempo.txt" | cut -d= -f2)/$P"
            AJ=$(grep '^carga_ajena_max=' "$D/tiempo.txt" | cut -d= -f2)
            RES=$(grep '^estado=' "$D/tiempo.txt" | tail -1 | cut -d= -f2)
        fi
        case "$ST" in
            PENDING|RUNNING|CONFIGURING|COMPLETING) NCURSO=$((NCURSO + 1)) ;;
            *) if [ "$RES" = "OK" ]; then NOK=$((NOK + 1)); else NMAL=$((NMAL + 1)); [ "$RES" = "-" ] && RES="$ST"; fi ;;
        esac
        printf "%-3s %-3s %-3s %-9s %-11s %-12s %8s %-7s %-11s  %s\n" "$I" "r$R" "$P" "$ID" "$ST" "$NODO" "$MIN" "$HIL" "${AJ:--}" "$RES"
    done < logs/jobs_A2.txt
    echo
    echo "Terminadas bien: $NOK   con problema: $NMAL   en cola o corriendo: $NCURSO"
    if [ "$NCURSO" -eq 0 ]; then
        echo "Todas terminaron. Siguiente paso:  python3 analisis/parteA2_resumen.py"
    else
        echo "La cadena corre una a la vez: las que dicen PENDING esperan a la anterior (es normal)."
    fi
    [ "$NMAL" -eq 0 ] || echo "Hay corridas con problema: mira logs/ex_A2_p<p>_r<r>_<Job>.out y mandame una captura."
fi

[ -s logs/jobs_A2d.txt ] || [ -s logs/jobs_A2.txt ] || echo "Todavia no hay corridas de la A.2."
