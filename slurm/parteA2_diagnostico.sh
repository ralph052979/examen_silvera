#!/bin/bash
# =============================================================================
#  parteA2_diagnostico.sh - Tarea A.2, paso 1: como pedir p nucleos FISICOS
#                           para p hilos OpenMP en Granado.
#
#  En el intento anterior, --cpus-per-task=4 dio solo 2 nucleos fisicos
#  (Slurm cuenta los hilos de hyperthreading como CPUs). Este diagnostico
#  prueba 4 formas de pedir los recursos con corridas CORTAS del caso base
#  (NELM = 4, 1-4 min cada una) y mide en que nucleos fisicos corren los hilos:
#
#    A : --ntasks=1 --cpus-per-task=p --hint=nomultithread     (la forma directa)
#    B : --ntasks=1 --cpus-per-task=p --threads-per-core=1
#    C : --ntasks=p --cpus-per-task=1 --hint=nomultithread     (asignacion igual a la A.1)
#        y un solo proceso: srun --ntasks=1 --cpus-per-task=p
#    D : --ntasks=1 --cpus-per-task=2p --hint=multithread      (2 hilos de hardware por
#        nucleo: cuesta el doble y solo cabe hasta p = 12 por el limite de 24 CPUs)
#
#  Corridas: A, B, C y D con p = 4; B y C con p = 24 (6 jobs, ~1 h-nucleo en total).
#
#  Uso (desde la carpeta del repositorio):   bash slurm/parteA2_diagnostico.sh
#  Despues:                                  python3 analisis/parteA2_diagnostico.py
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

PARTICION="${PARTICION:-phys_fast}"
REG=logs/jobs_A2d.txt
NODO_MALO=granado-compute-18
falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }

case "$PARTICION" in phys_fast|phys_general) ;; *) falla "particion no permitida por el examen: $PARTICION" ;; esac

echo "=== Comprobaciones previas ==="
for f in INCAR KPOINTS POSCAR; do
    [ -s "entradas/caso_base/$f" ] || falla "falta entradas/caso_base/$f"
done
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
grep -q "hilos_muestras" slurm/vasp_job.sh || falla "slurm/vasp_job.sh es una version anterior (corre bash slurm/parteA2_preparar.sh)"
[ -d .git ] || falla "no hay repositorio git"
if [ -n "$(git status --porcelain -- entradas slurm analisis README.md)" ]; then
    git status --short -- entradas slurm analisis README.md
    falla "hay cambios sin guardar en Git (corre primero: bash slurm/parteA2_preparar.sh)"
fi
mkdir -p logs corridas
if [ -s "$REG" ] && [ "${FORZAR:-0}" != "1" ]; then
    cat "$REG"
    falla "el diagnostico ya se envio. Resultado: python3 analisis/parteA2_diagnostico.py"
fi
echo "  [OK] caso base, POTCAR, librerias y Git en orden"

EXCL=""
if command -v sinfo >/dev/null 2>&1 && [ -n "$(sinfo -h -N -n "$NODO_MALO" -o %N 2>/dev/null)" ]; then
    EXCL="--exclude=$NODO_MALO"
fi
HASH=$(git rev-parse --short HEAD)

echo
echo "Se enviaran 6 corridas CORTAS de diagnostico (NELM = 4): A4 B4 C4 D4 B24 C24"
echo "  particion $PARTICION, commit $HASH ${EXCL:+(se excluye $NODO_MALO)}"
echo "  costo total estimado: ~1 h-nucleo; duracion: unos 10-15 min mas la espera en cola"
echo "  Es normal que alguna variante falle o sea rechazada: eso es lo que se esta midiendo."
read -r -p "Escribe s y Enter para enviar: " RESP
case "$RESP" in s|S) ;; *) echo "No se envio nada."; exit 0 ;; esac

for item in A:4 B:4 C:4 D:4 B:24 C:24; do
    V="${item%%:*}"; p="${item##*:}"
    case "$V" in
        A) OPTS="--ntasks=1 --cpus-per-task=$p --hint=nomultithread";  EXTRA="" ;;
        B) OPTS="--ntasks=1 --cpus-per-task=$p --threads-per-core=1";  EXTRA="" ;;
        C) OPTS="--ntasks=$p --cpus-per-task=1 --hint=nomultithread";  EXTRA=",PASO_NTASKS=1,PASO_CPT=$p" ;;
        D) OPTS="--ntasks=1 --cpus-per-task=$((2 * p)) --hint=multithread"; EXTRA="" ;;
    esac
    OUT=$(sbatch --parsable --job-name="ex_A2d_${V}${p}" --partition="$PARTICION" --nodes=1 $OPTS \
                 --mem=8G --time=00:20:00 $EXCL \
                 --export=ALL,CASO=caso_base,GIT_HASH="$HASH",NELM_DIAG=4,HILOS="$p"$EXTRA \
                 slurm/vasp_job.sh 2>&1)
    RC=$?
    ID=$(echo "$OUT" | tail -1)
    ID="${ID%%;*}"
    if [ "$RC" -eq 0 ] && echo "$ID" | grep -qE '^[0-9]+$'; then
        echo "$V $p $ID $(date '+%Y-%m-%d %H:%M')" >> "$REG"
        printf "  [OK] %s%-2s -> job %s\n" "$V" "$p" "$ID"
    else
        MSG=$(echo "$OUT" | tr '\n' ' ' | sed 's/  */ /g')
        echo "$V $p RECHAZADO $MSG" >> "$REG"
        printf "  [RECHAZADO POR SLURM] %s%-2s: %s\n" "$V" "$p" "$MSG"
    fi
done
echo
echo "Listo. En unos 15 min:  bash slurm/parteA2_estado.sh   y luego   python3 analisis/parteA2_diagnostico.py"
