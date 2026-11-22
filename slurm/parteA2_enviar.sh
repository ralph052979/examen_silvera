#!/bin/bash
# =============================================================================
#  parteA2_enviar.sh - Tarea A.2: escalamiento fuerte con OpenMP (memoria compartida).
#
#  1 proceso MPI con p hilos OpenMP, p = 2, 4, 8, 12, 16, 24 (hasta p_max), 3
#  repeticiones, en UNA sola cadena (--dependency=afterany). t(1) se toma de las
#  3 corridas con p = 1 de la A.1 (1 proceso y 1 hilo: es exactamente el mismo calculo).
#  La forma de pedir los recursos, el entorno de hilos (OPENBLAS_NUM_THREADS,
#  OMP_WAIT_POLICY) y el limite de tiempo de cada p salen del diagnostico
#  (datos/parteA2_config.txt, escrito por analisis/parteA2_hilos.py).
#
#  Uso (desde la carpeta del repositorio):
#      bash slurm/parteA2_enviar.sh               # envia el barrido
#      bash slurm/parteA2_enviar.sh --faltantes   # solo las corridas que faltan o fallaron
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

PARTICION="${PARTICION:-phys_fast}"
REG=logs/jobs_A2.txt
EXCLUIR=datos/parteA2_excluidas.txt
VARF=datos/parteA2_config.txt
MAX_PRIMERA="${MAX_PRIMERA:-420}"     # segundos; solo se cambia en pruebas
COSTO_MAX=40
NODO_MALO=granado-compute-18
ORDEN_R1="2 4 8 12 16 24"
ORDEN_R2="24 16 12 8 4 2"
ORDEN_R3="8 16 2 24 4 12"

tiempo_p() { grep "^limite_p$1=" "$VARF" | cut -d= -f2; }
falla()    { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }

case "$PARTICION" in phys_fast|phys_general) ;; *) falla "particion no permitida por el examen: $PARTICION" ;; esac
MODO=todas
case "${1:-}" in
    "") ;;
    --faltantes) MODO=faltantes ;;
    *) falla "opcion desconocida: $1 (usa --faltantes o nada)" ;;
esac

echo "=== Comprobaciones previas ==="
[ -s "$VARF" ] || falla "falta $VARF: primero la ronda 2 del diagnostico (bash slurm/parteA2_hilos.sh y python3 analisis/parteA2_hilos.py)"
git ls-files --error-unmatch "$VARF" >/dev/null 2>&1 || falla "$VARF no esta en Git: corre  python3 analisis/parteA2_hilos.py --guardar"
git diff --quiet HEAD -- "$VARF" || falla "$VARF cambio despues del commit: corre  python3 analisis/parteA2_hilos.py --guardar"
VAR=$(grep '^variante=' "$VARF" | cut -d= -f2)
PMAX=$(grep '^p_max=' "$VARF" | cut -d= -f2)
ENT=$(grep '^entorno=' "$VARF" | cut -d= -f2)
BLAS=$(grep '^blas_hilos=' "$VARF" | cut -d= -f2)
ESPERA=$(grep '^espera_omp=' "$VARF" | cut -d= -f2)
COSTO=$(grep '^costo_estimado_h=' "$VARF" | cut -d= -f2)
MODO_PASO=$(grep '^paso_modo=' "$VARF" | cut -d= -f2)
[ "$VAR" = "C" ] || falla "variante desconocida en $VARF: '$VAR' (vuelve a correr python3 analisis/parteA2_hilos.py --guardar)"
[ "$MODO_PASO" = "auto" ] || falla "paso_modo invalido en $VARF"
echo "$PMAX" | grep -qE '^[0-9]+$' || falla "p_max invalido en $VARF"
for p in 2 4 8 12 16 24; do
    [ "$p" -le "$PMAX" ] || continue
    tiempo_p "$p" | grep -qE '^[0-9]{2}:[0-9]{2}:[0-9]{2}$' || falla "falta limite_p$p en $VARF"
done
awk -v c="${COSTO:-999}" -v m="$COSTO_MAX" 'BEGIN { exit !(c + 0 <= m) }' || falla "el barrido costaria ${COSTO:-?} h-nucleo (mas de $COSTO_MAX): mandame una captura antes de enviarlo"
VARS_ENT=",MAX_PRIMERA_ITER=$MAX_PRIMERA,PASO_MODO=auto"
[ -n "$BLAS" ]       && VARS_ENT="$VARS_ENT,BLAS_HILOS=$BLAS"
[ -n "$ESPERA" ]     && VARS_ENT="$VARS_ENT,ESPERA_OMP=$ESPERA"
for f in INCAR KPOINTS POSCAR; do
    [ -s "entradas/caso_base/$f" ] || falla "falta entradas/caso_base/$f"
done
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
grep -q "PASO_MODO" slurm/vasp_job.sh || falla "slurm/vasp_job.sh es una version anterior (descomprime el paquete nuevo y corre bash slurm/parteA2_preparar.sh)"
if [ -n "$(git status --porcelain -- entradas slurm analisis README.md)" ]; then
    git status --short -- entradas slurm analisis README.md
    falla "hay cambios sin guardar en Git (corre primero: bash slurm/parteA2_preparar.sh)"
fi
mkdir -p logs corridas
if [ "$MODO" = "todas" ] && [ -s "$REG" ]; then
    falla "la A.2 ya se envio (ver $REG). Estado: bash slurm/parteA2_estado.sh  |  reenviar lo que falte: bash slurm/parteA2_enviar.sh --faltantes"
fi
if [ "$MODO" = "faltantes" ] && [ ! -s "$REG" ]; then
    falla "todavia no se ha enviado la A.2: usa  bash slurm/parteA2_enviar.sh"
fi
echo "  [OK] variante $VAR (p hasta $PMAX), entorno $ENT, caso base, POTCAR, librerias y Git en orden"

# ---- plan
PLAN=""
for r in 1 2 3; do
    eval "ORD=\$ORDEN_R$r"
    for p in $ORD; do
        [ "$p" -le "$PMAX" ] || continue
        if [ "$MODO" = "faltantes" ]; then
            hecho=0
            while read -r P R ID F H; do
                [ "$P" = "$p" ] && [ "$R" = "$r" ] || continue
                ST=$(squeue -h -j "$ID" -o %T 2>/dev/null | head -1)
                case "$ST" in PENDING|RUNNING|CONFIGURING|COMPLETING) hecho=1; continue ;; esac
                [ -f "$EXCLUIR" ] && grep -qE "^$ID([[:space:]]|$)" "$EXCLUIR" && continue
                E=$(grep '^estado=' "corridas/ex_A2_p${p}_r${r}_${ID}/tiempo.txt" 2>/dev/null | tail -1 | cut -d= -f2)
                [ "$E" = "OK" ] && hecho=1
            done < "$REG"
            [ "$hecho" -eq 1 ] && continue
        fi
        PLAN="$PLAN $r:$p"
    done
done
N=$(echo $PLAN | wc -w)
[ "$N" -gt 0 ] || { echo; echo "No falta ninguna corrida."; exit 0; }

EXCL=""
if command -v sinfo >/dev/null 2>&1 && [ -n "$(sinfo -h -N -n "$NODO_MALO" -o %N 2>/dev/null)" ]; then
    EXCL="--exclude=$NODO_MALO"
fi
HASH=$(git rev-parse --short HEAD)

echo
echo "Se enviaran $N corridas en una sola cadena (r = repeticion, p = hilos OpenMP de 1 proceso):"
echo "  $(echo $PLAN | sed 's/\([0-9]*\):\([0-9]*\)/r\1:p\2/g')"
echo "  nucleos como en la A.1 (--ntasks=p --hint=nomultithread) y paso de srun elegido por su mascara real"
echo "  entorno $ENT (OPENBLAS_NUM_THREADS=${BLAS:-sin fijar}, OMP_WAIT_POLICY=${ESPERA:-sin fijar})"
echo "  particion $PARTICION, commit $HASH ${EXCL:+(se excluye $NODO_MALO)}"
echo "  t(1) se toma de las 3 corridas con p = 1 de la A.1 (mismo calculo: 1 proceso, 1 hilo)"
if [ "$MODO" = "todas" ]; then
    echo "  costo estimado (proyeccion de la ronda 2): ${COSTO} h-nucleo"
fi
echo "  limites de tiempo: $(for p in 2 4 8 12 16 24; do [ "$p" -le "$PMAX" ] && printf "p%s=%s " "$p" "$(tiempo_p "$p")"; done)"
read -r -p "Escribe s y Enter para enviar: " RESP
case "$RESP" in s|S) ;; *) echo "No se envio nada."; exit 0 ;; esac

PREV=""
if [ "$MODO" = "faltantes" ]; then
    while read -r P R ID F H; do
        ST=$(squeue -h -j "$ID" -o %T 2>/dev/null | head -1)
        case "$ST" in PENDING|RUNNING|CONFIGURING|COMPLETING) PREV="$ID" ;; esac
    done < "$REG"
fi
for item in $PLAN; do
    r="${item%%:*}"; p="${item##*:}"
    OPTS="--ntasks=$p --cpus-per-task=1 --hint=nomultithread"
    DEP=""
    [ -n "$PREV" ] && DEP="--dependency=afterany:$PREV"
    ID=$(sbatch --parsable --job-name="ex_A2_p${p}_r${r}" --partition="$PARTICION" --nodes=1 $OPTS \
                --mem=8G --time="$(tiempo_p "$p")" $EXCL $DEP \
                --export=ALL,CASO=caso_base,GIT_HASH="$HASH",HILOS="$p"$VARS_ENT slurm/vasp_job.sh)
    [ -n "$ID" ] || falla "sbatch no devolvio un numero de job para r$r p=$p (las anteriores si quedaron enviadas; revisa con bash slurm/parteA2_estado.sh)"
    ID="${ID%%;*}"
    echo "$p $r $ID $(date '+%Y-%m-%d %H:%M')" >> "$REG"
    printf "  [OK] r%s  p = %-2s -> job %s\n" "$r" "$p" "$ID"
    PREV="$ID"
done
echo
echo "Listo. Seguimiento:  bash slurm/parteA2_estado.sh"
