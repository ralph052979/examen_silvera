#!/bin/bash
# =============================================================================
#  parteA2_enviar.sh - Tarea A.2: escalamiento fuerte con OpenMP (memoria compartida).
#
#  1 proceso MPI con p hilos OpenMP, p = 2, 4, 8, 12, 16, 24 (hasta p_max), 3
#  repeticiones, en UNA sola cadena (--dependency=afterany). t(1) se toma de las
#  3 corridas con p = 1 de la A.1 (1 proceso y 1 hilo: es exactamente el mismo calculo).
#  La forma de pedir los recursos es la que eligio el diagnostico
#  (datos/parteA2_variante.txt).
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
VARF=datos/parteA2_variante.txt
NODO_MALO=granado-compute-18
ORDEN_R1="2 4 8 12 16 24"
ORDEN_R2="24 16 12 8 4 2"
ORDEN_R3="8 16 2 24 4 12"

tiempo_p() { case "$1" in 2) echo 01:00:00 ;; 4) echo 00:45:00 ;; *) echo 00:30:00 ;; esac; }
falla()    { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }

case "$PARTICION" in phys_fast|phys_general) ;; *) falla "particion no permitida por el examen: $PARTICION" ;; esac
MODO=todas
case "${1:-}" in
    "") ;;
    --faltantes) MODO=faltantes ;;
    *) falla "opcion desconocida: $1 (usa --faltantes o nada)" ;;
esac

echo "=== Comprobaciones previas ==="
[ -s "$VARF" ] || falla "falta $VARF: primero el diagnostico (bash slurm/parteA2_diagnostico.sh y python3 analisis/parteA2_diagnostico.py)"
git ls-files --error-unmatch "$VARF" >/dev/null 2>&1 || falla "$VARF no esta en Git: corre  python3 analisis/parteA2_diagnostico.py --guardar"
git diff --quiet HEAD -- "$VARF" || falla "$VARF cambio despues del commit: corre  python3 analisis/parteA2_diagnostico.py --guardar"
VAR=$(grep '^variante=' "$VARF" | cut -d= -f2)
PMAX=$(grep '^p_max=' "$VARF" | cut -d= -f2)
case "$VAR" in B|C|D) ;; *) falla "variante desconocida en $VARF: '$VAR'" ;; esac
echo "$PMAX" | grep -qE '^[0-9]+$' || falla "p_max invalido en $VARF"
for f in INCAR KPOINTS POSCAR; do
    [ -s "entradas/caso_base/$f" ] || falla "falta entradas/caso_base/$f"
done
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
grep -q "hilos_muestras" slurm/vasp_job.sh || falla "slurm/vasp_job.sh es una version anterior (corre bash slurm/parteA2_preparar.sh)"
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
echo "  [OK] variante $VAR (p hasta $PMAX), caso base, POTCAR, librerias y Git en orden"

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
echo "  variante $VAR, particion $PARTICION, commit $HASH ${EXCL:+(se excluye $NODO_MALO)}"
echo "  t(1) se toma de las 3 corridas con p = 1 de la A.1 (mismo calculo: 1 proceso, 1 hilo)"
if [ "$MODO" = "todas" ]; then
    if [ "$VAR" = "D" ]; then
        echo "  costo estimado: unas 40 h-nucleo (la variante D cobra el doble); duracion: unas 3 h mas la cola"
    else
        echo "  costo estimado: unas 20-25 h-nucleo (se cobra el tiempo real); duracion: unas 3 h mas la cola"
    fi
fi
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
    case "$VAR" in
        C) OPTS="--ntasks=$p --cpus-per-task=1 --hint=nomultithread";  EXTRA=",PASO_NTASKS=1,PASO_CPT=$p" ;;
        B) OPTS="--ntasks=1 --cpus-per-task=$p --threads-per-core=1";  EXTRA="" ;;
        D) OPTS="--ntasks=1 --cpus-per-task=$((2 * p)) --hint=multithread"; EXTRA="" ;;
    esac
    DEP=""
    [ -n "$PREV" ] && DEP="--dependency=afterany:$PREV"
    ID=$(sbatch --parsable --job-name="ex_A2_p${p}_r${r}" --partition="$PARTICION" --nodes=1 $OPTS \
                --mem=8G --time="$(tiempo_p "$p")" $EXCL $DEP \
                --export=ALL,CASO=caso_base,GIT_HASH="$HASH",HILOS="$p"$EXTRA slurm/vasp_job.sh)
    [ -n "$ID" ] || falla "sbatch no devolvio un numero de job para r$r p=$p (las anteriores si quedaron enviadas; revisa con bash slurm/parteA2_estado.sh)"
    ID="${ID%%;*}"
    echo "$p $r $ID $(date '+%Y-%m-%d %H:%M')" >> "$REG"
    printf "  [OK] r%s  p = %-2s -> job %s\n" "$r" "$p" "$ID"
    PREV="$ID"
done
echo
echo "Listo. Seguimiento:  bash slurm/parteA2_estado.sh"
