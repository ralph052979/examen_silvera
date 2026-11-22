#!/bin/bash
# =============================================================================
#  enviar.sh - envia los candidatos de calibracion de la Tarea 0, cada uno en
#              1 nucleo fisico, SCF completo.
#
#  Uso (desde la carpeta del repositorio):
#      bash slurm/enviar.sh            # todos los candidatos de entradas/candidatos.txt
#      bash slurm/enviar.sh K8         # solo los que se nombran
#  Opcional:  PARTICION=phys_general bash slurm/enviar.sh
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

PARTICION="${PARTICION:-phys_fast}"
TIEMPO="${TIEMPO:-01:20:00}"
MEM="${MEM:-4G}"
REG=logs/jobs_t0.txt
NODO_MALO=granado-compute-18     # mataba los jobs en 1-2 s el 2026-11-18

falla() { echo "  [FALLA] $1"; echo; echo "No se envio nada mas. Mandame una captura."; exit 1; }

case "$PARTICION" in phys_fast|phys_general) ;; *) falla "particion no permitida por el examen: $PARTICION" ;; esac

if [ $# -gt 0 ]; then
    LISTA="$*"
else
    LISTA=$(awk '$1 !~ /^#/ && NF >= 4 {printf "%s ", $1}' entradas/candidatos.txt)
fi

echo "=== Comprobaciones previas ==="
[ -s entradas/POTCAR ]       || falla "falta entradas/POTCAR (corre primero: bash slurm/preparar.sh)"
[ -s lib/libfftw3_omp.so.3 ] || falla "falta lib/libfftw3_omp.so.3 (corre primero: bash slurm/preparar.sh)"
[ -s lib/libfftw3.so.3 ]     || falla "falta lib/libfftw3.so.3 (corre primero: bash slurm/preparar.sh)"
for K in $LISTA; do
    for f in INCAR KPOINTS POSCAR; do
        [ -s "entradas/candidatos/$K/$f" ] || falla "falta entradas/candidatos/$K/$f"
    done
done
[ -d .git ] || falla "no hay repositorio git (corre primero: bash slurm/preparar.sh)"
if [ -n "$(git status --porcelain -- entradas slurm analisis README.md)" ]; then
    git status --short -- entradas slurm analisis README.md
    falla "hay cambios sin guardar en Git. Las entradas deben estar en un commit antes de enviar (reproducibilidad)."
fi
mkdir -p logs corridas
if [ "${FORZAR:-0}" != "1" ] && [ -s "$REG" ]; then
    for K in $LISTA; do
        if grep -qE "^$K " "$REG"; then
            cat "$REG"
            falla "$K ya se envio (ver logs/jobs_t0.txt). Revisa con: bash slurm/estado.sh"
        fi
    done
fi
echo "  [OK] entradas, POTCAR, librerias y Git en orden"

EXCL=""
if command -v sinfo >/dev/null 2>&1 && [ -n "$(sinfo -h -N -n "$NODO_MALO" -o %N 2>/dev/null)" ]; then
    EXCL="--exclude=$NODO_MALO"
fi
HASH=$(git rev-parse --short HEAD)

echo
echo "Se enviaran: $LISTA"
echo "  particion $PARTICION, 1 nucleo fisico por job, --mem=$MEM, limite $TIEMPO ${EXCL:+(se excluye $NODO_MALO)}"
echo "  commit $HASH"
read -r -p "Escribe s y Enter para enviar: " R
case "$R" in s|S) ;; *) echo "No se envio nada."; exit 0 ;; esac

for K in $LISTA; do
    ID=$(sbatch --parsable --job-name="ex_T0b_$K" --partition="$PARTICION" \
                --nodes=1 --ntasks=1 --cpus-per-task=1 --hint=nomultithread \
                --mem="$MEM" --time="$TIEMPO" $EXCL \
                --export=ALL,CASO="$K",GIT_HASH="$HASH" slurm/vasp_job.sh)
    [ -n "$ID" ] || falla "sbatch no devolvio un numero de job para $K"
    ID="${ID%%;*}"
    echo "$K $ID $(date '+%Y-%m-%d %H:%M')" >> "$REG"
    echo "  [OK] $K -> job $ID"
done
echo
echo "Listo. Seguimiento:  bash slurm/estado.sh"
