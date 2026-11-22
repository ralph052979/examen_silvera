#!/bin/bash
# =============================================================================
#  parteA2_hilos.sh - Tarea A.2, ronda 2 del diagnostico de OpenMP: por que los
#                     hilos de VASP no rinden y con que entorno se corrige.
#
#  Hallazgo de la ronda 1: con p CPUs para p hilos, VASP no termino 4 iteraciones
#  en 20 min; con 2p CPUs (D4) termino, pero mas lento que con 1 hilo, y ps mostro
#  p + 1 hilos trabajando. Hipotesis: OpenBLAS (sin OpenMP) lanza su propio hilo,
#  que compite con los de OpenMP; como libgomp espera girando, quitarle el nucleo a
#  un solo hilo frena a todos.
#
#  Intento 1 (2026-11-22, logs/jobs_A2h.txt): con la variante B el paso de E1 p = 4
#  recibio 2 nucleos fisicos para 4 hilos y se detuvo al arrancar; las demas igual.
#  Intento 2 (este script, logs/jobs_A2h2.txt): el job pide los nucleos igual que la
#  A.1 (--ntasks=p --hint=nomultithread: p nucleos fisicos verificados, cobra p CPUs)
#  y vasp_job.sh (PASO_MODO=auto) prueba antes de VASP varias formas de lanzar 1
#  proceso con p hilos y usa la primera cuya mascara real tiene p nucleos fisicos.
#
#  Diseno 2 x 2 con p = 4:
#      E0 : entorno de la ronda 1 (control; se guarda la pila de los hilos con gdb)
#      E1 : OPENBLAS_NUM_THREADS=1                       (OpenBLAS sin hilos propios)
#      E2 : OMP_WAIT_POLICY=PASSIVE                      (OpenMP espera durmiendo)
#      E3 : OPENBLAS_NUM_THREADS=1 y OMP_WAIT_POLICY=PASSIVE
#  y con p = 24: E1 y E3, cada una solo si paso con p = 4 (si no, se omite en
#  segundos). Todo en UNA cadena (sin competencia entre jobs propios).
#  Cada corrida se detiene sola si a los 7 min no termino ninguna iteracion, o si la
#  mascara no tiene p nucleos fisicos (en segundos).
#
#  Uso (desde la carpeta del repositorio):   bash slurm/parteA2_hilos.sh
#  Despues:                                  python3 analisis/parteA2_hilos.py
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

PARTICION="${PARTICION:-phys_fast}"
REG=logs/jobs_A2h2.txt
NODO_MALO=granado-compute-18
MAX_PRIMERA="${MAX_PRIMERA:-420}"     # segundos; solo se cambia en pruebas
falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }

case "$PARTICION" in phys_fast|phys_general) ;; *) falla "particion no permitida por el examen: $PARTICION" ;; esac

echo "=== Comprobaciones previas ==="
for f in INCAR KPOINTS POSCAR; do
    [ -s "entradas/caso_base/$f" ] || falla "falta entradas/caso_base/$f"
done
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
grep -q "PASO_MODO" slurm/vasp_job.sh || falla "slurm/vasp_job.sh es una version anterior (descomprime el paquete nuevo y corre bash slurm/parteA2_preparar.sh)"
[ -d .git ] || falla "no hay repositorio git"
if [ -n "$(git status --porcelain -- entradas slurm analisis README.md)" ]; then
    git status --short -- entradas slurm analisis README.md
    falla "hay cambios sin guardar en Git (corre primero: bash slurm/parteA2_preparar.sh)"
fi
mkdir -p logs corridas
if [ -s "$REG" ] && [ "${FORZAR:-0}" != "1" ]; then
    cat "$REG"
    falla "el intento 2 de la ronda 2 ya se envio. Estado: bash slurm/parteA2_estado.sh  |  resultado: python3 analisis/parteA2_hilos.py"
fi
[ -s "$REG" ] && mv "$REG" "$REG.$(date +%Y%m%d%H%M%S)"
if command -v squeue >/dev/null 2>&1 && [ -s logs/jobs_A2h.txt ]; then
    for J in $(awk '$3 ~ /^[0-9]+$/ {print $3}' logs/jobs_A2h.txt); do
        [ -n "$(squeue -h -j "$J" -o %i 2>/dev/null)" ] && falla "el job $J del intento 1 sigue en cola o corriendo: espera a que termine (bash slurm/parteA2_estado.sh)"
    done
fi
echo "  [OK] caso base, POTCAR, librerias y Git en orden; nada del intento 1 sigue en cola"

EXCL=""
if command -v sinfo >/dev/null 2>&1 && [ -n "$(sinfo -h -N -n "$NODO_MALO" -o %N 2>/dev/null)" ]; then
    EXCL="--exclude=$NODO_MALO"
fi
HASH=$(git rev-parse --short HEAD)

echo
echo "Se enviaran 6 corridas CORTAS (NELM = 4) en una sola cadena (intento 2):"
echo "  p = 4 : E1 (OpenBLAS 1 hilo)  E3 (E1 + espera pasiva)  E2 (espera pasiva)  E0 (control)"
echo "  p = 24: E1 y E3, cada una solo si paso con p = 4 (si no, se omite en segundos)"
echo "  nucleos pedidos como en la A.1 (--ntasks=p --hint=nomultithread); el paso de srun se elige"
echo "  en cada corrida por su mascara real; particion $PARTICION, commit $HASH ${EXCL:+(se excluye $NODO_MALO)}"
echo "  costo estimado: 3-6 h-nucleo; duracion: 30-45 min mas la espera en cola"
echo "  Es normal que E0 (y quiza E2) falle: es el control que reproduce el problema."
read -r -p "Escribe s y Enter para enviar: " RESP
case "$RESP" in s|S) ;; *) echo "No se envio nada."; exit 0 ;; esac

PREV=""
declare -A ID
enviar() {   # enviar <etiqueta> <p> <variables extra>
    local E="$1" p="$2" VARS="$3" DEP="" OUT RC J
    [ -n "$PREV" ] && DEP="--dependency=afterany:$PREV"
    OUT=$(sbatch --parsable --job-name="ex_A2h_${E}_p${p}" --partition="$PARTICION" --nodes=1 \
                 --ntasks="$p" --cpus-per-task=1 --hint=nomultithread \
                 --mem=8G --time=00:15:00 $EXCL $DEP \
                 --export=ALL,CASO=caso_base,GIT_HASH="$HASH",NELM_DIAG=4,HILOS="$p",PASO_MODO=auto,MAX_PRIMERA_ITER=$MAX_PRIMERA$VARS \
                 slurm/vasp_job.sh 2>&1)
    RC=$?
    J=$(echo "$OUT" | tail -1); J="${J%%;*}"
    if [ "$RC" -ne 0 ] || ! echo "$J" | grep -qE '^[0-9]+$'; then
        echo "$E $p RECHAZADO $(echo "$OUT" | tr '\n' ' ')" >> "$REG"
        falla "Slurm rechazo ${E} p=$p: $(echo "$OUT" | tr '\n' ' ') (las anteriores si quedaron enviadas)"
    fi
    echo "$E $p $J $(date '+%Y-%m-%d %H:%M')" >> "$REG"
    printf "  [OK] %s p = %-2s -> job %s\n" "$E" "$p" "$J"
    ID["$E$p"]="$J"
    PREV="$J"
}
enviar E1 4  ",BLAS_HILOS=1"
enviar E3 4  ",BLAS_HILOS=1,ESPERA_OMP=PASSIVE"
enviar E2 4  ",ESPERA_OMP=PASSIVE"
enviar E0 4  ",PILA_HILOS=12"
enviar E1 24 ",BLAS_HILOS=1,REQUIERE_OK=${ID[E14]}"
enviar E3 24 ",BLAS_HILOS=1,ESPERA_OMP=PASSIVE,REQUIERE_OK=${ID[E34]}"
echo
echo "Listo. Seguimiento:  bash slurm/parteA2_estado.sh    Resultado:  python3 analisis/parteA2_hilos.py"
