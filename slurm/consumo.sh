#!/bin/bash
# =============================================================================
#  consumo.sh - horas-nucleo consumidas (registro de consumo que pide el examen).
#  Uso:  bash slurm/consumo.sh [FECHA_INICIO]     (por defecto 2026-11-16)
#  Guarda la salida completa de sacct en datos/sacct_consumo.txt
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
DESDE="${1:-2026-11-16}"
mkdir -p datos

sacct -u "$USER" -S "$DESDE" -X \
      --format=JobID,JobName%22,Partition,AllocCPUS,Elapsed,CPUTimeRAW,MaxRSS,State \
      > datos/sacct_consumo.txt 2>&1 \
    || { echo "sacct fallo:"; cat datos/sacct_consumo.txt; exit 1; }

sacct -u "$USER" -S "$DESDE" -X -n -P \
      --format=JobID,JobName,Partition,AllocCPUS,Elapsed,CPUTimeRAW,State 2>/dev/null |
awk -F'|' -v desde="$DESDE" '
    $3 == "phys_fast" || $3 == "phys_general" {
        h = $6 / 3600.0; tot += h; n++
        if ($2 ~ /^ex_/)     { ex += h; nex++ }
        if ($2 ~ /^ex_T0b_/) { t0 += h; nt0++ }
    }
    END {
        printf "Jobs en phys_fast y phys_general desde %s: %d  ->  %.2f h-nucleo\n", desde, n, tot
        printf "  del examen (nombre ex_*):                    %d  ->  %.2f h-nucleo\n", nex, ex
        printf "  de esta Tarea 0 (nombre ex_T0b_*):           %d  ->  %.2f h-nucleo\n", nt0, t0
        printf "Presupuesto restante: %.1f de 250 h-nucleo (contando todo phys_fast y phys_general)\n", 250 - tot
    }'
echo "Detalle completo guardado en datos/sacct_consumo.txt"
