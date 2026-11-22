#!/bin/bash
# =============================================================================
#  parteC_enviar.sh - Tarea C.2: precision vs costo con ENCUT.
#
#  ENCUT = 250, 300, 350, 400, 450, 500, 550, 600, 700, 800 eV (referencia: 800) con la
#  malla 3x3x1 (entradas/encut/, generadas por parteBC_preparar.sh), SCF completo, en
#  1 nucleo (--ntasks=1 --hint=nomultithread, como la p = 1 de la A.1).
#  3 repeticiones; cada repeticion es UN job que corre la serie completa, una corrida
#  tras otra (slurm/parteC2_serie.sh), en orden distinto en cada repeticion. Asi solo
#  se agregan 3 jobs a la cola (limite de 50 por usuario).
#  Los 3 jobs van en cadena y empiezan cuando termina la ultima corrida de la Parte B
#  que siga en cola (--dependency=afterany): nunca corren a la vez que la Parte B.
#
#  Uso (desde la carpeta del repositorio):
#      bash slurm/parteC_enviar.sh               # envia las 3 series
#      bash slurm/parteC_enviar.sh --faltantes   # solo las corridas que faltan o fallaron
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

PARTICION="${PARTICION:-phys_fast}"
REG=logs/jobs_C2.txt
EXCLUIR=datos/parteC_excluidas.txt
MAX_PRIMERA="${MAX_PRIMERA_C2:-300}"     # segundos; solo se cambia en pruebas
LIMITE=03:00:00
MAX_COLA=50
NODO_MALO=granado-compute-18
ORDEN_R1="250 300 350 400 450 500 550 600 700 800"
ORDEN_R2="800 700 600 550 500 450 400 350 300 250"
ORDEN_R3="500 250 700 400 800 300 600 350 550 450"

falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }

case "$PARTICION" in phys_fast|phys_general) ;; *) falla "particion no permitida por el examen: $PARTICION" ;; esac
MODO=todas
case "${1:-}" in
    "") ;;
    --faltantes) MODO=faltantes ;;
    *) falla "opcion desconocida: $1 (usa --faltantes o nada)" ;;
esac

echo "=== Comprobaciones previas ==="
for E in $ORDEN_R1; do
    for f in INCAR KPOINTS POSCAR; do
        [ -s "entradas/encut/E$E/$f" ] || falla "falta entradas/encut/E$E/$f (corre primero: bash slurm/parteBC_preparar.sh)"
    done
    [ "$(grep -c "^ENCUT  = $E\$" "entradas/encut/E$E/INCAR")" -eq 1 ] || falla "entradas/encut/E$E/INCAR no tiene ENCUT = $E"
done
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
[ -s slurm/parteC2_serie.sh ] && grep -q 'SUFIJO:-' slurm/vasp_job.sh \
    || falla "faltan scripts nuevos (descomprime el paquete nuevo y corre bash slurm/parteBC_preparar.sh)"
if [ -n "$(git status --porcelain -- entradas slurm analisis README.md)" ]; then
    git status --short -- entradas slurm analisis README.md
    falla "hay cambios sin guardar en Git (corre primero: bash slurm/parteBC_preparar.sh)"
fi
git ls-files --error-unmatch entradas/encut/E800/INCAR >/dev/null 2>&1 \
    || falla "las entradas de la C.2 no estan en Git (corre primero: bash slurm/parteBC_preparar.sh)"
mkdir -p logs corridas
if [ "$MODO" = "todas" ] && [ -s "$REG" ]; then
    falla "la C.2 ya se envio (ver $REG). Estado: bash slurm/parteBC_estado.sh  |  reenviar lo que falte: bash slurm/parteC_enviar.sh --faltantes"
fi
if [ "$MODO" = "faltantes" ] && [ ! -s "$REG" ]; then
    falla "todavia no se ha enviado la C.2: usa  bash slurm/parteC_enviar.sh"
fi
EN_COLA=0; PREV=""
if command -v squeue >/dev/null 2>&1; then
    EN_COLA=$(squeue -h -u "$(id -un)" -o %i 2>/dev/null | wc -l)
    # la cadena nueva sigue a la ultima corrida de la Parte B o de la C.2 que siga en cola
    PREV=$(squeue -h -u "$(id -un)" -o "%i %j" 2>/dev/null | awk '$2 ~ /^ex_(B1|B2|C2)_/ {print $1}' | sort -n | tail -1)
fi
echo "  [OK] entradas, POTCAR, librerias y Git en orden"

en_cola_id() {
    case "$(squeue -h -j "$1" -o %T 2>/dev/null | head -1)" in PENDING|RUNNING|CONFIGURING|COMPLETING) return 0 ;; esac
    return 1
}
corrida_ok() {   # corrida_ok <r> <E>: 0 si alguna corrida de esa repeticion y ENCUT termino OK
    local r="$1" E="$2" R ID RESTO
    while read -r R ID RESTO; do
        [ "$R" = "$r" ] || continue
        [ -f "$EXCLUIR" ] && grep -qE "^($ID|${ID}_E$E)([[:space:]]|$)" "$EXCLUIR" && continue
        [ "$(grep '^estado=' "corridas/ex_C2_r${r}_${ID}_E$E/tiempo.txt" 2>/dev/null | tail -1 | cut -d= -f2)" = "OK" ] && return 0
    done < "$REG"
    return 1
}
PLAN=""
for r in 1 2 3; do
    eval "ORD=\$ORDEN_R$r"
    LISTA=""
    if [ "$MODO" = "faltantes" ]; then
        ocupada=0
        while read -r R ID RESTO; do
            [ "$R" = "$r" ] && en_cola_id "$ID" && ocupada=1
        done < "$REG"
        [ "$ocupada" = 1 ] && continue
        for E in $ORD; do corrida_ok "$r" "$E" || LISTA="$LISTA:$E"; done
        LISTA="${LISTA#:}"
    else
        LISTA=$(echo $ORD | tr ' ' ':')
    fi
    [ -n "$LISTA" ] && PLAN="$PLAN $r/$LISTA"
done
N=$(echo $PLAN | wc -w)
[ "$N" -gt 0 ] || { echo; echo "No falta ninguna corrida (o las que faltan siguen en cola)."; exit 0; }
[ $((EN_COLA + N)) -le "$MAX_COLA" ] \
    || falla "ya tienes $EN_COLA jobs en cola y se agregarian $N (maximo $MAX_COLA por usuario): espera a que terminen algunos"

EXCL=""
if command -v sinfo >/dev/null 2>&1 && [ -n "$(sinfo -h -N -n "$NODO_MALO" -o %N 2>/dev/null)" ]; then
    EXCL="--exclude=$NODO_MALO"
fi
HASH=$(git rev-parse --short HEAD)

echo
echo "Se enviaran $N jobs de 1 nucleo en cadena (cada uno corre una serie de ENCUT, malla 3x3x1, SCF completo):"
for item in $PLAN; do echo "  r${item%%/*}: ENCUT = $(echo "${item#*/}" | tr ':' ' ') eV"; done
if [ -n "$PREV" ]; then
    echo "  empiezan cuando termine el job $PREV (la ultima corrida en cola de las Partes B/C)"
fi
echo "  limite $LIMITE por serie; particion $PARTICION, commit $HASH ${EXCL:+(se excluye $NODO_MALO)}"
[ "$MODO" = "todas" ] && echo "  costo estimado: unas 3 h-nucleo en total; cada serie dura alrededor de 1 h"
read -r -p "Escribe s y Enter para enviar: " RESP
case "$RESP" in s|S) ;; *) echo "No se envio nada."; exit 0 ;; esac

for item in $PLAN; do
    r="${item%%/*}"; LISTA="${item#*/}"
    DEP=""
    [ -n "$PREV" ] && DEP="--dependency=afterany:$PREV"
    ID=$(sbatch --parsable --job-name="ex_C2_r$r" --partition="$PARTICION" --nodes=1 \
                --ntasks=1 --cpus-per-task=1 --hint=nomultithread \
                --mem=4G --time="$LIMITE" $EXCL $DEP \
                --export=ALL,ENCUTS="$LISTA",GIT_HASH="$HASH",BLAS_HILOS=1,MAX_PRIMERA_ITER="$MAX_PRIMERA" \
                slurm/parteC2_serie.sh)
    [ -n "$ID" ] || falla "sbatch no devolvio un numero de job para la serie r$r (las anteriores si quedaron enviadas; revisa con bash slurm/parteBC_estado.sh)"
    ID="${ID%%;*}"
    echo "$r $ID $LISTA $(date '+%Y-%m-%d %H:%M')" >> "$REG"
    printf "  [OK] C.2 r%s -> job %s\n" "$r" "$ID"
    PREV="$ID"
done
echo
echo "Listo. Seguimiento:  bash slurm/parteBC_estado.sh"
