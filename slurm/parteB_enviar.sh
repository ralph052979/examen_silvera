#!/bin/bash
# =============================================================================
#  parteB_enviar.sh - Parte B: escalamiento debil (B.1) e hibrido MPI x OpenMP (B.2).
#
#  B.1  p procesos MPI (p = 1, 2, 4, 8, 12, 16, 24), KPAR = p, 3 k-points irreducibles
#       por proceso (entradas/debil/P<p>/, generadas por parteBC_preparar.sh), SCF
#       completo. Nucleos y lanzamiento como en la A.1 (--ntasks=p --hint=nomultithread,
#       srun --cpu-bind=cores). Ew(p) = t(1)/t(p).
#  B.2  24 nucleos fijos (--ntasks=24 --hint=nomultithread, como la A.1 y la A.2) y el
#       paso de srun con N procesos MPI x T hilos OpenMP:
#       (24,1) (12,2) (8,3) (6,4) (4,6) (2,12), 3 repeticiones, y (1,24) 1 vez (sus otras
#       3 repeticiones son las p = 24 de la A.2, con la misma carga y el mismo entorno).
#       Carga fija de 4 iteraciones SCF del caso base (NELM = 4), como la A.2: t(1) es la
#       p = 1 de la A.2. El paso se elige por su mascara real (PASO_MODO=auto): cada
#       proceso con T nucleos fisicos propios; OMP_PLACES=cores, OMP_PROC_BIND=close.
#  En todas: OPENBLAS_NUM_THREADS=1 y memoria maxima (VmHWM) de cada proceso cada 3 s.
#
#  Todo va en UNA cadena (--dependency=afterany: nunca corren dos a la vez), intercalando
#  B.1 y B.2 por repeticion y con orden distinto en cada repeticion.
#
#  Uso (desde la carpeta del repositorio):
#      bash slurm/parteB_enviar.sh               # envia la Parte B
#      bash slurm/parteB_enviar.sh --faltantes   # solo las corridas que faltan o fallaron
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

PARTICION="${PARTICION:-phys_fast}"
REG1=logs/jobs_B1.txt
REG2=logs/jobs_B2.txt
EXCLUIR=datos/parteB_excluidas.txt
MAX_PRIMERA_B1="${MAX_PRIMERA_B1:-240}"   # segundos; solo se cambian en pruebas
MAX_PRIMERA_B2="${MAX_PRIMERA_B2:-420}"
MON_B="${MON_B:-3}"                       # s entre muestras de hilos y memoria
NELM_FIJO=4
LIMITE_B1=00:30:00
LIMITE_B2=00:20:00
MAX_COLA=50                               # jobs en cola por usuario (regla del examen)
NODO_MALO=granado-compute-18
B1_R1="1 2 4 8 12 16 24"
B1_R2="24 16 12 8 4 2 1"
B1_R3="8 16 1 24 4 12 2"
B2_R1="24x1 12x2 8x3 6x4 4x6 2x12"
B2_R2="2x12 4x6 1x24 6x4 8x3 12x2 24x1"
B2_R3="6x4 24x1 2x12 8x3 4x6 12x2"
VARS_B1=",HILOS=1,BLAS_HILOS=1,MAX_PRIMERA_ITER=$MAX_PRIMERA_B1,MON_INTERVALO=$MON_B"
VARS_B2=",PASO_MODO=auto,NELM_DIAG=$NELM_FIJO,BLAS_HILOS=1,MAX_PRIMERA_ITER=$MAX_PRIMERA_B2,MON_INTERVALO=$MON_B"

falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }

case "$PARTICION" in phys_fast|phys_general) ;; *) falla "particion no permitida por el examen: $PARTICION" ;; esac
MODO=todas
case "${1:-}" in
    "") ;;
    --faltantes) MODO=faltantes ;;
    *) falla "opcion desconocida: $1 (usa --faltantes o nada)" ;;
esac

echo "=== Comprobaciones previas ==="
for f in INCAR KPOINTS POSCAR; do
    [ -s "entradas/caso_base/$f" ] || falla "falta entradas/caso_base/$f"
done
for p in $B1_R1; do
    for f in INCAR KPOINTS POSCAR; do
        [ -s "entradas/debil/P$p/$f" ] || falla "falta entradas/debil/P$p/$f (corre primero: bash slurm/parteBC_preparar.sh)"
    done
    [ "$(grep -c "^KPAR   = $p\$" "entradas/debil/P$p/INCAR")" -eq 1 ] || falla "entradas/debil/P$p/INCAR no tiene KPAR = $p"
done
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
grep -q "mem_total_mb" slurm/vasp_job.sh && grep -q 'NTASKS="${PASO_NTASKS:-1}"' slurm/vasp_job.sh \
    || falla "slurm/vasp_job.sh es una version anterior (descomprime el paquete nuevo y corre bash slurm/parteBC_preparar.sh)"
if [ -n "$(git status --porcelain -- entradas slurm analisis README.md)" ]; then
    git status --short -- entradas slurm analisis README.md
    falla "hay cambios sin guardar en Git (corre primero: bash slurm/parteBC_preparar.sh)"
fi
git ls-files --error-unmatch entradas/debil/P24/INCAR >/dev/null 2>&1 \
    || falla "las entradas de la B.1 no estan en Git (corre primero: bash slurm/parteBC_preparar.sh)"
git ls-files --error-unmatch datos/parteA2_resumen.csv >/dev/null 2>&1 \
    || falla "falta guardar la A.2: corre primero  python3 analisis/parteA2_resumen.py --guardar"
[ -s logs/jobs_A2.txt ] || falla "no encuentro logs/jobs_A2.txt (la B.2 usa la p = 1 y la p = 24 de la A.2)"
mkdir -p logs corridas
if [ "$MODO" = "todas" ] && { [ -s "$REG1" ] || [ -s "$REG2" ]; }; then
    falla "la Parte B ya se envio (ver $REG1 y $REG2). Estado: bash slurm/parteBC_estado.sh  |  reenviar lo que falte: bash slurm/parteB_enviar.sh --faltantes"
fi
if [ "$MODO" = "faltantes" ] && [ ! -s "$REG1" ] && [ ! -s "$REG2" ]; then
    falla "todavia no se ha enviado la Parte B: usa  bash slurm/parteB_enviar.sh"
fi
EN_COLA=0
if command -v squeue >/dev/null 2>&1; then
    [ -z "$(squeue -h -u "$(id -un)" -o %j 2>/dev/null | grep '^ex_A2')" ] \
        || falla "todavia hay corridas de la A.2 en cola: espera a que terminen (bash slurm/parteA2_estado.sh)"
    EN_COLA=$(squeue -h -u "$(id -un)" -o %i 2>/dev/null | wc -l)
fi
echo "  [OK] entradas, POTCAR, librerias, Git y A.2 guardada en orden"

# ---- plan: B1:p:r y B2:n:t:r, intercalados por repeticion
hecha() {   # hecha <registro> <nombre de job> <clave en el registro>: 1 si ya esta OK o en cola
    local reg="$1" nom="$2" clave="$3" h=0 L ID ST E
    [ -s "$reg" ] || { echo 0; return; }
    while read -r L; do
        case "$L" in "$clave "*) ;; *) continue ;; esac
        ID=$(echo "$L" | awk -v n="$(echo "$clave" | wc -w)" '{print $(n + 1)}')
        ST=$(squeue -h -j "$ID" -o %T 2>/dev/null | head -1)
        case "$ST" in PENDING|RUNNING|CONFIGURING|COMPLETING) h=1; continue ;; esac
        [ -f "$EXCLUIR" ] && grep -qE "^$ID([[:space:]]|$)" "$EXCLUIR" && continue
        E=$(grep '^estado=' "corridas/${nom}_${ID}/tiempo.txt" 2>/dev/null | tail -1 | cut -d= -f2)
        [ "$E" = "OK" ] && h=1
    done < "$reg"
    echo "$h"
}
PLAN=""
for r in 1 2 3; do
    eval "O1=\$B1_R$r; O2=\$B2_R$r"
    for p in $O1; do
        [ "$MODO" = "faltantes" ] && [ "$(hecha "$REG1" "ex_B1_p${p}_r${r}" "$p $r")" = 1 ] && continue
        PLAN="$PLAN B1:$p:$r"
    done
    for nt in $O2; do
        n="${nt%x*}"; t="${nt#*x}"; rr="$r"
        [ "$nt" = "1x24" ] && rr=1          # (1,24): una sola corrida nueva
        [ "$MODO" = "faltantes" ] && [ "$(hecha "$REG2" "ex_B2_n${n}x${t}_r${rr}" "$n $t $rr")" = 1 ] && continue
        PLAN="$PLAN B2:$n:$t:$rr"
    done
done
N=$(echo $PLAN | wc -w)
[ "$N" -gt 0 ] || { echo; echo "No falta ninguna corrida."; exit 0; }
N1=$(echo $PLAN | tr ' ' '\n' | grep -c '^B1:')
N2=$(echo $PLAN | tr ' ' '\n' | grep -c '^B2:')
[ $((EN_COLA + N)) -le "$MAX_COLA" ] \
    || falla "ya tienes $EN_COLA jobs en cola y se agregarian $N (maximo $MAX_COLA por usuario): espera a que terminen algunos"

EXCL=""
if command -v sinfo >/dev/null 2>&1 && [ -n "$(sinfo -h -N -n "$NODO_MALO" -o %N 2>/dev/null)" ]; then
    EXCL="--exclude=$NODO_MALO"
fi
HASH=$(git rev-parse --short HEAD)

echo
echo "Se enviaran $N corridas en una sola cadena ($N1 de la B.1 y $N2 de la B.2):"
echo "  B.1: p procesos MPI, KPAR = p, 3p k-points (3 por proceso), SCF completo, limite $LIMITE_B1"
echo "  B.2: 24 nucleos, N procesos x T hilos, 4 iteraciones SCF del caso base, limite $LIMITE_B2"
echo "  orden: $(echo $PLAN | sed -e 's/B1:\([0-9]*\):\([0-9]*\)/B1-p\1-r\2/g' -e 's/B2:\([0-9]*\):\([0-9]*\):\([0-9]*\)/B2-\1x\2-r\3/g')"
echo "  OPENBLAS_NUM_THREADS=1, memoria de cada proceso cada 3 s; particion $PARTICION, commit $HASH ${EXCL:+(se excluye $NODO_MALO)}"
if [ "$MODO" = "todas" ]; then
    echo "  costo estimado: unas 40 h-nucleo (B.1 unas 17, B.2 unas 22); duracion: unas 5 h mas la espera en cola"
fi
read -r -p "Escribe s y Enter para enviar: " RESP
case "$RESP" in s|S) ;; *) echo "No se envio nada."; exit 0 ;; esac

PREV=""
if [ "$MODO" = "faltantes" ]; then     # la nueva cadena sigue a la ultima corrida que siga en cola
    for ID in $(awk '{print $3}' "$REG1" 2>/dev/null) $(awk '{print $4}' "$REG2" 2>/dev/null); do
        ST=$(squeue -h -j "$ID" -o %T 2>/dev/null | head -1)
        case "$ST" in
            PENDING|RUNNING|CONFIGURING|COMPLETING)
                if [ -z "$PREV" ] || [ "$ID" -gt "$PREV" ]; then PREV="$ID"; fi ;;
        esac
    done
fi
for item in $PLAN; do
    DEP=""
    [ -n "$PREV" ] && DEP="--dependency=afterany:$PREV"
    case "$item" in
        B1:*)
            IFS=: read -r _ p r <<< "$item"
            NOM="ex_B1_p${p}_r${r}"
            ID=$(sbatch --parsable --job-name="$NOM" --partition="$PARTICION" --nodes=1 \
                        --ntasks="$p" --cpus-per-task=1 --hint=nomultithread \
                        --mem=$((4 + p))G --time="$LIMITE_B1" $EXCL $DEP \
                        --export=ALL,CASO=debil/P$p,GIT_HASH="$HASH"$VARS_B1 slurm/vasp_job.sh)
            [ -n "$ID" ] || falla "sbatch no devolvio un numero de job para $NOM (las anteriores si quedaron enviadas; revisa con bash slurm/parteBC_estado.sh)"
            ID="${ID%%;*}"
            echo "$p $r $ID $(date '+%Y-%m-%d %H:%M')" >> "$REG1"
            printf "  [OK] B.1 r%s  p = %-2s      -> job %s\n" "$r" "$p" "$ID"
            ;;
        B2:*)
            IFS=: read -r _ n t r <<< "$item"
            NOM="ex_B2_n${n}x${t}_r${r}"
            ID=$(sbatch --parsable --job-name="$NOM" --partition="$PARTICION" --nodes=1 \
                        --ntasks=24 --cpus-per-task=1 --hint=nomultithread \
                        --mem=24G --time="$LIMITE_B2" $EXCL $DEP \
                        --export=ALL,CASO=caso_base,GIT_HASH="$HASH",HILOS="$t",PASO_NTASKS="$n"$VARS_B2 slurm/vasp_job.sh)
            [ -n "$ID" ] || falla "sbatch no devolvio un numero de job para $NOM (las anteriores si quedaron enviadas; revisa con bash slurm/parteBC_estado.sh)"
            ID="${ID%%;*}"
            echo "$n $t $r $ID $(date '+%Y-%m-%d %H:%M')" >> "$REG2"
            printf "  [OK] B.2 r%s  %2s x %-2s    -> job %s\n" "$r" "$n" "$t" "$ID"
            ;;
    esac
    PREV="$ID"
done
echo
echo "Listo. Seguimiento:  bash slurm/parteBC_estado.sh"
