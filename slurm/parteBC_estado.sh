#!/bin/bash
# =============================================================================
#  parteBC_estado.sh - estado de las corridas de la Parte B y de la Tarea C.2
#                      (no gasta horas-nucleo).
#  Uso:  bash slurm/parteBC_estado.sh
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

estado_job() {   # imprime "ESTADO MINUTOS"
    local L
    L=$(squeue -h -j "$1" -o "%T %M" 2>/dev/null | head -1)
    if [ -n "$L" ]; then echo "$L"; else echo "$(sacct -n -X -j "$1" -o State%20 2>/dev/null | head -1 | tr -d ' ') -"; fi
}
NC=0; NOK=0; NMAL=0; NTOT=0
fila() {   # fila <etiqueta> <jobid> <directorio>
    local ET="$1" ID="$2" D="$3" ST MIN NODO="-" IT="-" MEM="-" RES="-" W
    read -r ST MIN <<< "$(estado_job "$ID")"
    [ -n "$ST" ] || ST="?"
    NTOT=$((NTOT + 1))
    if [ ! -s "$D/tiempo.txt" ] || ! grep -q '^estado=' "$D/tiempo.txt"; then
        case "$ST" in PENDING|RUNNING|CONFIGURING|COMPLETING) NC=$((NC + 1)) ;; *) NMAL=$((NMAL + 1)) ;; esac
    fi
    IT=$( { wc -l < "$D/progreso.txt"; } 2>/dev/null || echo -)
    if [ -s "$D/tiempo.txt" ]; then
        RES=$(grep '^estado=' "$D/tiempo.txt" | tail -1 | cut -d= -f2)
        NODO=$(grep '^nodo=' "$D/entorno.txt" 2>/dev/null | cut -d= -f2 | cut -d. -f1 | sed 's/granado-//')
        W=$(grep '^wall_s=' "$D/tiempo.txt" | cut -d= -f2)
        [ -n "$W" ] && MIN=$(awk -v w="$W" 'BEGIN{printf "%.1f", w/60}')
        M=$(grep '^mem_total_mb=' "$D/tiempo.txt" | cut -d= -f2)
        [ -n "$M" ] && [ "$M" != "0" ] && MEM=$(awk -v m="$M" 'BEGIN{printf "%.2f", m/1024}')
        case "$RES" in OK) NOK=$((NOK + 1)) ;; "") ;; *) NMAL=$((NMAL + 1)) ;; esac
    fi
    printf "%-14s %-9s %-11s %-12s %8s %5s %8s  %s\n" "$ET" "$ID" "$ST" "${NODO:--}" "${MIN:--}" "${IT:--}" "$MEM" "$(echo "$RES" | cut -c1-60)"
}
cab() { printf "%-14s %-9s %-11s %-12s %8s %5s %8s  %s\n" Corrida Job Slurm Nodo Minutos Iter "Mem(GB)" Resultado; }

if [ -s logs/jobs_B1.txt ]; then
    echo "=== B.1 escalamiento debil (p procesos MPI, KPAR = p, 3p k-points, SCF completo) ==="
    cab
    while read -r P R ID F H; do
        [ -n "$ID" ] || continue
        fila "p$P r$R" "$ID" "corridas/ex_B1_p${P}_r${R}_${ID}"
    done < logs/jobs_B1.txt
    echo
fi
if [ -s logs/jobs_B2.txt ]; then
    echo "=== B.2 hibrido (24 nucleos: N procesos x T hilos, 4 iteraciones SCF) ==="
    cab
    while read -r N T R ID F H; do
        [ -n "$ID" ] || continue
        fila "${N}x${T} r$R" "$ID" "corridas/ex_B2_n${N}x${T}_r${R}_${ID}"
    done < logs/jobs_B2.txt
    echo
fi
NB=$NTOT; NCB=$NC
if [ -s logs/jobs_C2.txt ]; then
    echo "=== C.2 serie de ENCUT (malla 3x3x1, 1 nucleo, SCF completo; un job por repeticion) ==="
    cab
    while read -r R ID LISTA F H; do
        [ -n "$ID" ] || continue
        read -r STJ MINJ <<< "$(estado_job "$ID")"
        for E in $(echo "$LISTA" | tr ':' ' '); do
            D="corridas/ex_C2_r${R}_${ID}_E$E"
            if [ -d "$D" ]; then
                fila "E$E r$R" "$ID" "$D"
            else
                NTOT=$((NTOT + 1))
                case "$STJ" in PENDING|RUNNING|CONFIGURING|COMPLETING) NC=$((NC + 1)); ST="${STJ}" ;; *) ST="NO_CORRIO"; NMAL=$((NMAL + 1)) ;; esac
                printf "%-14s %-9s %-11s %-12s %8s %5s %8s  %s\n" "E$E r$R" "$ID" "$ST" - - - - -
            fi
        done
    done < logs/jobs_C2.txt
    echo
fi
if [ "$NTOT" -eq 0 ]; then
    echo "Todavia no se ha enviado la Parte B: bash slurm/parteB_enviar.sh"
    exit 0
fi
echo "Resumen: $NTOT corridas | $NOK OK | $NMAL con otro resultado | $NC en cola o corriendo"
echo "(Iter = iteraciones SCF terminadas; Mem = suma de la memoria maxima de los procesos de VASP)"
if [ "$NB" -gt 0 ] && [ "$NCB" -eq 0 ]; then
    echo "La Parte B termino. Siguiente paso:  python3 analisis/parteB_resumen.py"
fi
if [ ! -s logs/jobs_C2.txt ]; then
    echo "La C.2 no se ha enviado:  bash slurm/parteC_enviar.sh"
elif [ "$NC" -eq 0 ]; then
    echo "Todas terminaron (B y C.2). Mandame una captura de este resumen."
fi
