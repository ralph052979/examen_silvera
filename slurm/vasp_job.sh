#!/bin/bash
#SBATCH --partition=phys_fast
#SBATCH --nodes=1
#SBATCH --hint=nomultithread
#SBATCH --output=logs/%x_%j.out
#
# =============================================================================
#  vasp_job.sh - ejecuta UN calculo VASP en su propia carpeta
#                corridas/<nombre_del_job>_<jobid>/
#
#  No se envia a mano: lo envia slurm/enviar.sh, que fija --job-name,
#  --ntasks, --cpus-per-task, --mem, --time y la variable CASO.
#
#  Variables de entrada (sbatch --export=ALL,...):
#    CASO      K1..K7     -> entradas/candidatos/<CASO>/  (INCAR KPOINTS POSCAR)
#              caso_base  -> entradas/caso_base/
#    GIT_HASH  (opcional) commit del repositorio al momento de enviar
#
#  Configuracion de Granado verificada el 2026-11-17 (ver README.md):
#    * --hint=nomultithread    : 1 proceso por nucleo FISICO; AllocCPUS = ntasks
#    * srun --mpi=pmix         : unico lanzador que respeta la asignacion de Slurm
#                                (mpirun ata todo al nucleo 0; srun solo falla con PMIx)
#    * lib/ en LD_LIBRARY_PATH : libfftw3_omp.so.3 no existe en los nodos de computo
#    * tiempo con date         : /usr/bin/time no existe en los nodos
#    * NO fijar PMIX_MCA_gds   : rompe MPI_Init. Los avisos "PMIX ERROR: BAD-PARAM"
#                                son inofensivos.
# =============================================================================

set -u
REPO="${SLURM_SUBMIT_DIR:?este script se envia con sbatch desde la carpeta del repositorio}"
cd "$REPO" || exit 1

terminar() {   # terminar <codigo> <estado>
    echo "estado=$2" >> "${RUN:-$REPO/logs}/tiempo.txt" 2>/dev/null
    echo "RESULTADO: $2"
    exit "$1"
}

# ------------------------------------------------------------- 1. entradas
CASO="${CASO:-}"
[ -n "$CASO" ] || { echo "RESULTADO: FALLA_falta_la_variable_CASO"; exit 1; }
if [ "$CASO" = "caso_base" ]; then
    SRC="$REPO/entradas/caso_base"
else
    SRC="$REPO/entradas/candidatos/$CASO"
fi
for f in INCAR KPOINTS POSCAR; do
    [ -s "$SRC/$f" ] || { echo "RESULTADO: FALLA_no_existe_$f (en $SRC)"; exit 1; }
done
[ -s "$REPO/entradas/POTCAR" ] || { echo "RESULTADO: FALLA_no_existe_entradas/POTCAR"; exit 1; }

RUN="$REPO/corridas/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir -p "$RUN" || exit 1
cp "$SRC/INCAR" "$SRC/KPOINTS" "$SRC/POSCAR" "$REPO/entradas/POTCAR" "$RUN/"

# ------------------------------------------------------------- 2. entorno
if ! type module >/dev/null 2>&1; then
    for f in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
             /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash; do
        [ -f "$f" ] && . "$f" && break
    done
fi
module purge
module load vasp/6.3.0-gnu9-openmp
export LD_LIBRARY_PATH="$REPO/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
ulimit -s unlimited 2>/dev/null || echo "aviso: no se pudo fijar ulimit -s unlimited"

cd "$RUN" || exit 1
VASP_BIN="$(command -v vasp_std)" || terminar 2 "FALLA_vasp_std_no_encontrado"

{
    echo "caso=$CASO"
    echo "job_id=$SLURM_JOB_ID"
    echo "job_name=$SLURM_JOB_NAME"
    echo "nodo=$(hostname)"
    echo "particion=${SLURM_JOB_PARTITION:-?}"
    echo "ntasks=${SLURM_NTASKS:-1}"
    echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-1}"
    echo "omp_num_threads=$OMP_NUM_THREADS"
    echo "git_hash=${GIT_HASH:-desconocido}"
    echo "vasp_bin=$VASP_BIN"
    echo "fecha_inicio=$(date '+%Y-%m-%d %H:%M:%S')"
} > entorno.txt
module list > modulos.txt 2>&1
ldd "$VASP_BIN" > ldd.txt 2>&1
env | grep -E '^(SLURM_|OMP_|LD_LIBRARY_PATH)' | sort > variables.txt
scontrol show job -d "$SLURM_JOB_ID" > slurm_job.txt 2>&1
md5sum INCAR KPOINTS POSCAR POTCAR > entradas.md5

if grep -q "not found" ldd.txt; then
    grep "not found" ldd.txt
    terminar 2 "FALLA_librerias_faltantes"
fi

# ------------------------------------------------------------- 3. carga del nodo
# Mide cuantos CPUs logicos estan ocupados (>50 %) en el nodo y en el socket del job.
# Sirve para saber si otros usuarios compiten por el mismo socket (ver README).
PRIMER_CPU=$(grep -o "CPU_IDs=[0-9]*" slurm_job.txt | head -1 | cut -d= -f2)
PKG=$(cat "/sys/devices/system/cpu/cpu${PRIMER_CPU:-0}/topology/physical_package_id" 2>/dev/null || echo 0)
SOCKET_CPUS=""
for d in /sys/devices/system/cpu/cpu[0-9]*; do
    n="${d##*cpu}"
    p=$(cat "$d/topology/physical_package_id" 2>/dev/null || echo 0)
    [ "$p" = "$PKG" ] && SOCKET_CPUS="$SOCKET_CPUS $n"
done
carga() {
    grep '^cpu[0-9]' /proc/stat > .stat1 2>/dev/null
    sleep 2
    grep '^cpu[0-9]' /proc/stat > .stat2 2>/dev/null
    awk -v sock="$SOCKET_CPUS" -v etq="$1" -v f="$(date '+%H:%M:%S')" -v la="$(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null)" '
        BEGIN { n = split(sock, s, " "); for (i = 1; i <= n; i++) en["cpu" s[i]] = 1 }
        FNR == NR { t = 0; for (i = 2; i <= 9; i++) t += $i; T[$1] = t; I[$1] = $5 + $6; next }
        { t = 0; for (i = 2; i <= 9; i++) t += $i; dt = t - T[$1]; di = ($5 + $6) - I[$1]
          if (dt > 0) { u = 1 - di / dt; tot++; if (u > 0.5) oc++
                        if ($1 in en) { ts++; if (u > 0.5) os++ } } }
        END { printf "%s %s cpus_ocupadas_nodo=%d/%d cpus_ocupadas_socket=%d/%d loadavg=%s\n",
                     etq, f, oc, tot, os, ts, la }' .stat1 .stat2
    rm -f .stat1 .stat2
}
carga inicio > carga.txt
( while sleep 300; do carga durante >> carga.txt; done ) 2>/dev/null &
MONITOR=$!

# ------------------------------------------------------------- 4. VASP
T0=$(date +%s.%N)
srun --mpi=pmix --cpu-bind=verbose,cores vasp_std > vasp.out 2> vasp.err
RC=$?
T1=$(date +%s.%N)
pkill -P "$MONITOR" 2>/dev/null
kill "$MONITOR" 2>/dev/null
wait "$MONITOR" 2>/dev/null
carga fin >> carga.txt
WALL=$(awk -v a="$T0" -v b="$T1" 'BEGIN{printf "%.3f", b-a}')

# ------------------------------------------------------------- 5. verificacion
NTASKS="${SLURM_NTASKS:-1}"
NMASK=$(grep -o "mask 0x[0-9a-fA-F]*" vasp.err 2>/dev/null | sort -u | wc -l)
FIN=$(grep -c "General timing and accounting" OUTCAR 2>/dev/null); FIN=${FIN:-0}
CONV=$(grep -c "aborting loop because EDIFF is reached" OUTCAR 2>/dev/null); CONV=${CONV:-0}
NSCF=$(grep -c "LOOP:" OUTCAR 2>/dev/null); NSCF=${NSCF:-0}

{
    echo "fecha_fin=$(date '+%Y-%m-%d %H:%M:%S')"
    echo "wall_s=$WALL"
    echo "rc_srun=$RC"
    echo "mascaras_distintas=$NMASK"
    echo "fin_normal_vasp=$FIN"
    echo "convergido=$CONV"
    echo "n_scf=$NSCF"
} > tiempo.txt

echo "---- afinidad (srun --cpu-bind=verbose) ----"
grep "mask" vasp.err
echo "---- carga del nodo ----"
cat carga.txt
echo "---- tiempo ----"
cat tiempo.txt

[ "$RC" -eq 0 ]            || { tail -20 vasp.out; tail -20 vasp.err; terminar 3 "FALLA_srun_rc_$RC"; }
[ "$FIN" -ge 1 ]           || { tail -20 vasp.out; terminar 3 "FALLA_VASP_no_termino_normalmente"; }
[ "$NSCF" -ge 1 ]          || terminar 3 "FALLA_sin_iteraciones_SCF"
[ "$NMASK" -eq "$NTASKS" ] || terminar 4 "FALLA_afinidad_${NMASK}_mascaras_para_${NTASKS}_procesos"
[ "$CONV" -ge 1 ]          || terminar 5 "FALLA_SCF_no_convergio"
terminar 0 "OK"
