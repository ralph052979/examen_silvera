#!/bin/bash
# =============================================================================
#  parteA1_enviar.sh - Tarea A.1: escalamiento fuerte con MPI.
#
#  21 corridas del caso base: p = 1, 2, 4, 8, 12, 16, 24 procesos MPI x 3
#  repeticiones, en UNA sola cadena (cada job espera a que termine el anterior,
#  --dependency=afterany). Asi ninguna corrida propia compite con otra por el
#  nodo. Orden: r1 ascendente, r2 descendente, r3 mezclado, para que la carga
#  del cluster a distintas horas no se confunda con el efecto de p.
#
#  Uso (desde la carpeta del repositorio):
#      bash slurm/parteA1_enviar.sh               # envia las 21 corridas
#      bash slurm/parteA1_enviar.sh --faltantes   # solo las que faltan o fallaron
#  Opcional:  PARTICION=phys_general bash slurm/parteA1_enviar.sh
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

PARTICION="${PARTICION:-phys_fast}"
REG=logs/jobs_A1.txt
EXCLUIR=datos/parteA1_excluidas.txt
NODO_MALO=granado-compute-18     # mataba los jobs en 1-2 s el 2026-11-18
ORDEN_R1="1 2 4 8 12 16 24"
ORDEN_R2="24 16 12 8 4 2 1"
ORDEN_R3="8 1 16 4 24 2 12"

tiempo_p() { case "$1" in 1) echo 01:40:00 ;; 2) echo 01:00:00 ;; 4) echo 00:40:00 ;; 8) echo 00:30:00 ;; *) echo 00:20:00 ;; esac; }
mem_p()    { echo "$((4 + $1 / 2))G"; }
falla()    { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }

case "$PARTICION" in phys_fast|phys_general) ;; *) falla "particion no permitida por el examen: $PARTICION" ;; esac
MODO=todas
case "${1:-}" in
    "") ;;
    --faltantes) MODO=faltantes ;;
    *) falla "opcion desconocida: $1 (usa --faltantes o nada)" ;;
esac

echo "=== Comprobaciones previas ==="
for f in INCAR KPOINTS POSCAR; do
    [ -s "entradas/caso_base/$f" ] || falla "falta entradas/caso_base/$f (la Tarea 0 no se fijo)"
done
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
grep -q "afinidad_ok" slurm/vasp_job.sh || falla "slurm/vasp_job.sh es la version anterior (corre bash slurm/parteA1_preparar.sh)"
[ -d .git ] || falla "no hay repositorio git"
if [ -n "$(git status --porcelain -- entradas slurm analisis README.md)" ]; then
    git status --short -- entradas slurm analisis README.md
    falla "hay cambios sin guardar en Git (corre primero: bash slurm/parteA1_preparar.sh)"
fi
mkdir -p logs corridas
if [ "$MODO" = "todas" ] && [ -s "$REG" ]; then
    falla "la A.1 ya se envio (ver $REG). Estado: bash slurm/parteA1_estado.sh  |  reenviar lo que falte: bash slurm/parteA1_enviar.sh --faltantes"
fi
if [ "$MODO" = "faltantes" ] && [ ! -s "$REG" ]; then
    falla "todavia no se ha enviado la A.1: usa  bash slurm/parteA1_enviar.sh"
fi
echo "  [OK] caso base, POTCAR, librerias y Git en orden"

# ---- plan: lista de r:p en el orden de ejecucion
PLAN=""
for r in 1 2 3; do
    eval "ORD=\$ORDEN_R$r"
    for p in $ORD; do
        if [ "$MODO" = "faltantes" ]; then
            hecho=0
            while read -r P R ID F H; do
                [ "$P" = "$p" ] && [ "$R" = "$r" ] || continue
                ST=$(squeue -h -j "$ID" -o %T 2>/dev/null | head -1)
                case "$ST" in PENDING|RUNNING|CONFIGURING|COMPLETING) hecho=1; continue ;; esac
                [ -f "$EXCLUIR" ] && grep -qE "^$ID([[:space:]]|$)" "$EXCLUIR" && continue
                E=$(grep '^estado=' "corridas/ex_A1_p${p}_r${r}_${ID}/tiempo.txt" 2>/dev/null | tail -1 | cut -d= -f2)
                [ "$E" = "OK" ] && hecho=1
            done < "$REG"
            [ "$hecho" -eq 1 ] && continue
        fi
        PLAN="$PLAN $r:$p"
    done
done
N=$(echo $PLAN | wc -w)
[ "$N" -gt 0 ] || { echo; echo "No falta ninguna corrida: las 21 terminaron bien (o siguen en cola)."; exit 0; }

EXCL=""
if command -v sinfo >/dev/null 2>&1 && [ -n "$(sinfo -h -N -n "$NODO_MALO" -o %N 2>/dev/null)" ]; then
    EXCL="--exclude=$NODO_MALO"
fi
HASH=$(git rev-parse --short HEAD)

echo
echo "Se enviaran $N corridas en una sola cadena, en este orden (r = repeticion, p = procesos MPI):"
echo "  $(echo $PLAN | sed 's/\([0-9]*\):\([0-9]*\)/r\1:p\2/g')"
echo "  particion $PARTICION, 1 nucleo fisico por proceso, NCORE = 1, KPAR = 1 ${EXCL:+(se excluye $NODO_MALO)}"
echo "  commit $HASH"
if [ "$MODO" = "todas" ]; then
    echo "  costo estimado: unas 14 h-nucleo (se cobra el tiempo real); duracion: unas 4 h mas la espera en cola"
fi
read -r -p "Escribe s y Enter para enviar: " RESP
case "$RESP" in s|S) ;; *) echo "No se envio nada."; exit 0 ;; esac

# en --faltantes, la cadena nueva espera a la ultima corrida que siga en cola
PREV=""
if [ "$MODO" = "faltantes" ]; then
    while read -r P R ID F H; do
        ST=$(squeue -h -j "$ID" -o %T 2>/dev/null | head -1)
        case "$ST" in PENDING|RUNNING|CONFIGURING|COMPLETING) PREV="$ID" ;; esac
    done < "$REG"
fi
for item in $PLAN; do
    r="${item%%:*}"; p="${item##*:}"
    DEP=""
    [ -n "$PREV" ] && DEP="--dependency=afterany:$PREV"
    ID=$(sbatch --parsable --job-name="ex_A1_p${p}_r${r}" --partition="$PARTICION" \
                --nodes=1 --ntasks="$p" --cpus-per-task=1 --hint=nomultithread \
                --mem="$(mem_p "$p")" --time="$(tiempo_p "$p")" $EXCL $DEP \
                --export=ALL,CASO=caso_base,GIT_HASH="$HASH" slurm/vasp_job.sh)
    [ -n "$ID" ] || falla "sbatch no devolvio un numero de job para r$r p=$p (las anteriores si quedaron enviadas; revisa con bash slurm/parteA1_estado.sh)"
    ID="${ID%%;*}"
    echo "$p $r $ID $(date '+%Y-%m-%d %H:%M')" >> "$REG"
    printf "  [OK] r%s  p = %-2s -> job %s\n" "$r" "$p" "$ID"
    PREV="$ID"
done
echo
echo "Listo. Seguimiento:  bash slurm/parteA1_estado.sh"
