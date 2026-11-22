#!/bin/bash
# =============================================================================
#  parteA2_enviar.sh - Tarea A.2: escalamiento fuerte con OpenMP (memoria compartida).
#
#  1 proceso MPI con p hilos OpenMP, p = 1, 2, 4, 8, 12, 16, 24, 3 repeticiones, en
#  UNA sola cadena (--dependency=afterany), con orden intercalado.
#
#  Carga fija de 4 iteraciones SCF (NELM = 4), no el SCF completo: el diagnostico
#  (ronda 2, datos/parteA2_hilos.md) mostro que con 24 hilos una iteracion tarda MAS
#  que con 1 hilo, asi que cada corrida completa con p hilos costaria casi p x 30 min
#  (unas 110 h-nucleo el barrido). Con 4 iteraciones el barrido cuesta unas 31 h-nucleo.
#  t(1) se mide con la misma carga (p = 1, tambien 3 veces). Para comparar con MPI se
#  usa el tiempo por iteracion (LOOP del OUTCAR), que la A.1 tambien tiene.
#
#  Entorno (ronda 2): OPENBLAS_NUM_THREADS=1 (E1; igual de rapido que sin fijarlo y
#  asegura que solo VASP lanza hilos) y la espera de libgomp por defecto (la espera
#  pasiva fue 50-70 % mas lenta). Nucleos pedidos como en la A.1
#  (--ntasks=p --hint=nomultithread) y paso de srun elegido por su mascara real
#  (PASO_MODO=auto, ver geometria_pasos.txt de cada corrida).
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
MAX_PRIMERA="${MAX_PRIMERA:-420}"     # segundos; solo se cambia en pruebas
NELM_FIJO=4
LIMITE=00:25:00
NODO_MALO=granado-compute-18
ORDEN_R1="1 2 4 8 12 16 24"
ORDEN_R2="24 16 12 8 4 2 1"
ORDEN_R3="8 16 1 24 4 12 2"
VARS_ENT=",NELM_DIAG=$NELM_FIJO,BLAS_HILOS=1,PASO_MODO=auto,MAX_PRIMERA_ITER=$MAX_PRIMERA"

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
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
grep -q "PASO_MODO" slurm/vasp_job.sh || falla "slurm/vasp_job.sh es una version anterior (descomprime el paquete nuevo y corre bash slurm/parteA2_preparar.sh)"
grep -q "carga fija de 4 iteraciones" analisis/parteA2_resumen.py || falla "analisis/parteA2_resumen.py es una version anterior (descomprime el paquete nuevo)"
if [ -n "$(git status --porcelain -- entradas slurm analisis README.md)" ]; then
    git status --short -- entradas slurm analisis README.md
    falla "hay cambios sin guardar en Git (corre primero: bash slurm/parteA2_preparar.sh)"
fi
git ls-files --error-unmatch datos/parteA2_hilos.md >/dev/null 2>&1 \
    || falla "falta guardar el diagnostico: corre primero  python3 analisis/parteA2_hilos.py --guardar"
mkdir -p logs corridas
if [ "$MODO" = "todas" ] && [ -s "$REG" ]; then
    falla "la A.2 ya se envio (ver $REG). Estado: bash slurm/parteA2_estado.sh  |  reenviar lo que falte: bash slurm/parteA2_enviar.sh --faltantes"
fi
if [ "$MODO" = "faltantes" ] && [ ! -s "$REG" ]; then
    falla "todavia no se ha enviado la A.2: usa  bash slurm/parteA2_enviar.sh"
fi
if command -v squeue >/dev/null 2>&1 && [ -n "$(squeue -h -u "$(id -un)" -n ex_A2h_E1_p4,ex_A2h_E3_p4,ex_A2h_E2_p4,ex_A2h_E0_p4,ex_A2h_E1_p24,ex_A2h_E3_p24 -o %i 2>/dev/null)" ]; then
    falla "todavia hay corridas del diagnostico en cola: espera a que terminen (bash slurm/parteA2_estado.sh)"
fi
echo "  [OK] caso base, POTCAR, librerias, Git y diagnostico guardado en orden"

# ---- plan
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
echo "  carga fija: 4 iteraciones SCF (NELM = $NELM_FIJO) del caso base; t(1) con la misma carga"
echo "  nucleos como en la A.1 (--ntasks=p --hint=nomultithread), paso de srun elegido por su mascara real"
echo "  OPENBLAS_NUM_THREADS=1, espera de OpenMP por defecto; limite $LIMITE por corrida"
echo "  particion $PARTICION, commit $HASH ${EXCL:+(se excluye $NODO_MALO)}"
if [ "$MODO" = "todas" ]; then
    echo "  costo estimado: unas 31 h-nucleo; duracion: unas 3.5 h mas la espera en cola"
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
    DEP=""
    [ -n "$PREV" ] && DEP="--dependency=afterany:$PREV"
    ID=$(sbatch --parsable --job-name="ex_A2_p${p}_r${r}" --partition="$PARTICION" --nodes=1 \
                --ntasks="$p" --cpus-per-task=1 --hint=nomultithread \
                --mem=8G --time="$LIMITE" $EXCL $DEP \
                --export=ALL,CASO=caso_base,GIT_HASH="$HASH",HILOS="$p"$VARS_ENT slurm/vasp_job.sh)
    [ -n "$ID" ] || falla "sbatch no devolvio un numero de job para r$r p=$p (las anteriores si quedaron enviadas; revisa con bash slurm/parteA2_estado.sh)"
    ID="${ID%%;*}"
    echo "$p $r $ID $(date '+%Y-%m-%d %H:%M')" >> "$REG"
    printf "  [OK] r%s  p = %-2s -> job %s\n" "$r" "$p" "$ID"
    PREV="$ID"
done
echo
echo "Listo. Seguimiento:  bash slurm/parteA2_estado.sh"
