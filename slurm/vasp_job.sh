#!/bin/bash
#SBATCH --partition=phys_fast
#SBATCH --nodes=1
#SBATCH --output=logs/%x_%j.out
#
# =============================================================================
#  vasp_job.sh - ejecuta UN calculo VASP en su propia carpeta
#                corridas/<nombre_del_job>_<jobid>/
#
#  No se envia a mano: lo envian slurm/enviar.sh (Tarea 0),
#  slurm/parteA1_enviar.sh (A.1), slurm/parteA2_diagnostico.sh y
#  slurm/parteA2_enviar.sh (A.2). Ellos fijan --job-name, --ntasks,
#  --cpus-per-task, --hint, --mem, --time y las variables de abajo.
#
#  Variables de entrada (sbatch --export=ALL,...):
#    CASO         K1..K7     -> entradas/candidatos/<CASO>/  (INCAR KPOINTS POSCAR)
#                 caso_base  -> entradas/caso_base/
#    GIT_HASH     (opcional) commit del repositorio al momento de enviar
#    PASO_NTASKS  (opcional) procesos del paso srun   (por defecto, los del job)
#    PASO_CPT     (opcional) CPUs por proceso del paso (por defecto, las del job)
#    HILOS        (opcional) hilos OpenMP por proceso  (por defecto PASO_CPT)
#    NELM_DIAG    (opcional) fija NELM: corrida corta de diagnostico (no exige convergencia)
#
#  Configuracion de Granado verificada (ver README.md):
#    * --hint=nomultithread    : 1 proceso por nucleo FISICO; AllocCPUS = ntasks
#                                (lo pasa cada script de envio, no esta aqui)
#    * srun --mpi=pmix         : unico lanzador que respeta la asignacion de Slurm
#    * lib/ en LD_LIBRARY_PATH : libfftw3_omp.so.3 no existe en los nodos de computo
#    * tiempo con date         : /usr/bin/time no existe en los nodos
#    * NO fijar PMIX_MCA_gds   : rompe MPI_Init. Los avisos "PMIX ERROR: BAD-PARAM"
#                                son inofensivos.
#
#  Verificaciones de cada corrida (quedan en tiempo.txt, afinidad.txt e hilos.txt):
#    * afinidad: cada proceso tiene exactamente HILOS nucleos FISICOS en su mascara
#      y ningun nucleo fisico se comparte entre procesos;
#    * hilos (si HILOS > 1): los hilos de VASP corren en PASO_NTASKS x HILOS nucleos
#      fisicos distintos (se mira con ps -L durante la corrida);
#    * VASP termino normalmente y el SCF convergio (salvo con NELM_DIAG);
#    * carga ajena: CPUs ocupados por otros usuarios en los sockets del job.
# =============================================================================

set -u
REPO="${SLURM_SUBMIT_DIR:?este script se envia con sbatch desde la carpeta del repositorio}"
cd "$REPO" || exit 1
SYSCPU="${SYSCPU:-/sys/devices/system/cpu}"    # solo se cambia en pruebas
PROCSTAT="${PROCSTAT:-/proc/stat}"              # solo se cambia en pruebas

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
if [ -n "${NELM_DIAG:-}" ]; then
    sed -i "s/^NELM .*/NELM   = ${NELM_DIAG}/" "$RUN/INCAR"
    [ "$(grep -c "^NELM   = ${NELM_DIAG}\$" "$RUN/INCAR")" -eq 1 ] || terminar 1 "FALLA_no_se_pudo_fijar_NELM"
fi

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
NTASKS="${PASO_NTASKS:-${SLURM_NTASKS:-1}}"
CPT="${PASO_CPT:-${SLURM_CPUS_PER_TASK:-1}}"
NHILOS="${HILOS:-$CPT}"
GEOM=""
[ -n "${PASO_NTASKS:-}" ] && GEOM="$GEOM --ntasks=$PASO_NTASKS"
[ -n "${PASO_CPT:-}" ]    && GEOM="$GEOM --cpus-per-task=$PASO_CPT"
export OMP_NUM_THREADS="$NHILOS"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
ulimit -s unlimited 2>/dev/null || echo "aviso: no se pudo fijar ulimit -s unlimited"

cd "$RUN" || exit 1
VASP_BIN="$(command -v vasp_std)" || terminar 2 "FALLA_vasp_std_no_encontrado"
CPU_MODELO=$(grep -m1 "model name" /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ *//')

{
    echo "caso=$CASO"
    echo "job_id=$SLURM_JOB_ID"
    echo "job_name=$SLURM_JOB_NAME"
    echo "nodo=$(hostname)"
    echo "cpu=${CPU_MODELO:-?}"
    echo "particion=${SLURM_JOB_PARTITION:-?}"
    echo "job_ntasks=${SLURM_NTASKS:-1}"
    echo "job_cpus_per_task=${SLURM_CPUS_PER_TASK:-1}"
    echo "paso_ntasks=$NTASKS"
    echo "paso_cpus_per_task=$CPT"
    echo "omp_num_threads=$OMP_NUM_THREADS"
    echo "nelm_diag=${NELM_DIAG:-no}"
    echo "git_hash=${GIT_HASH:-desconocido}"
    echo "vasp_bin=$VASP_BIN"
    echo "fecha_inicio=$(date '+%Y-%m-%d %H:%M:%S')"
} > entorno.txt
module list > modulos.txt 2>&1
ldd "$VASP_BIN" > ldd.txt 2>&1
env | grep -E '^(SLURM_|OMP_|LD_LIBRARY_PATH)' | sort > variables.txt
scontrol show job -d "$SLURM_JOB_ID" > slurm_job.txt 2>&1
ALLOC=$(grep -o "NumCPUs=[0-9]*" slurm_job.txt | head -1 | cut -d= -f2)
md5sum INCAR KPOINTS POSCAR POTCAR > entradas.md5

if grep -q "not found" ldd.txt; then
    grep "not found" ldd.txt
    terminar 2 "FALLA_librerias_faltantes"
fi

# ------------------------------------------------------------- 3. carga del nodo
# CPUs logicos ocupados (>50 %) en cada socket del nodo, medidos en 2 s.
MAPA=""
for d in "$SYSCPU"/cpu[0-9]*; do
    n="${d##*cpu}"
    p=$(cat "$d/topology/physical_package_id" 2>/dev/null || echo 0)
    MAPA="$MAPA $n:$p"
done
carga() {
    grep '^cpu[0-9]' "$PROCSTAT" > .stat1 2>/dev/null
    sleep 2
    grep '^cpu[0-9]' "$PROCSTAT" > .stat2 2>/dev/null
    awk -v mapa="$MAPA" -v etq="$1" -v f="$(date '+%H:%M:%S')" -v la="$(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null)" '
        BEGIN { n = split(mapa, m, " "); mx = -1
                for (i = 1; i <= n; i++) { split(m[i], a, ":"); pk["cpu" a[1]] = a[2]; tot[a[2]]++
                                           if (a[2] + 0 > mx) mx = a[2] + 0 } }
        FNR == NR { t = 0; for (i = 2; i <= 9; i++) t += $i; T[$1] = t; I[$1] = $5 + $6; next }
        { t = 0; for (i = 2; i <= 9; i++) t += $i; dt = t - T[$1]; di = ($5 + $6) - I[$1]
          if (dt > 0 && ($1 in pk)) { u = 1 - di / dt; if (u > 0.5) oc[pk[$1]]++ } }
        END { s = ""; for (k = 0; k <= mx; k++) if (k in tot) s = s sprintf(" s%d=%d/%d", k, oc[k] + 0, tot[k])
              printf "%s %s%s loadavg=%s\n", etq, f, s, la }' .stat1 .stat2
    rm -f .stat1 .stat2
}
carga inicio > carga.txt
( while sleep 60; do carga durante >> carga.txt; done ) 2>/dev/null &
MONITOR=$!

# ------------------------------------------------------------- 3b. monitor de hilos (OpenMP)
# Cada 10 s anota en que CPU logico esta cada hilo de los procesos vasp_std de ESTE job.
USUARIO="$(id -un)"
monitor_hilos() {
    sleep 15
    while :; do
        for pid in $(ps -u "$USUARIO" -o pid=,comm= 2>/dev/null | awk '$2 == "vasp_std" {print $1}'); do
            tr '\0' '\n' < "/proc/$pid/environ" 2>/dev/null | grep -qx "SLURM_JOB_ID=$SLURM_JOB_ID" || continue
            echo "muestra $(date '+%H:%M:%S') pid=$pid $(ps -L -o psr= -p "$pid" 2>/dev/null | tr -s ' \n' ' ')"
        done
        sleep 10
    done
}
: > hilos_crudo.txt
MONH=""
if [ "$NHILOS" -gt 1 ]; then
    monitor_hilos >> hilos_crudo.txt 2>/dev/null &
    MONH=$!
fi

# ------------------------------------------------------------- 4. VASP
T0=$(date +%s.%N)
srun --mpi=pmix $GEOM --cpu-bind=verbose,cores vasp_std > vasp.out 2> vasp.err
RC=$?
T1=$(date +%s.%N)
for m in $MONITOR $MONH; do
    pkill -P "$m" 2>/dev/null
    kill "$m" 2>/dev/null
    wait "$m" 2>/dev/null
done
carga fin >> carga.txt
WALL=$(awk -v a="$T0" -v b="$T1" 'BEGIN{printf "%.3f", b-a}')

# ------------------------------------------------------------- 5. afinidad por nucleo fisico
# srun --cpu-bind=verbose escribe en vasp.err una linea por proceso:
#   cpu-bind=MASK - <nodo>, task  3  3 [12345]: mask 0x8000000008 set
grep -E "task +[0-9]+ +[0-9]+ \[[0-9]+\]: mask 0x[0-9a-fA-F]+" vasp.err 2>/dev/null |
awk -v syscpu="$SYSCPU" -v cpt="$NHILOS" -v nt="$NTASKS" '
    function hex2cpus(h,    s, n, i, d, b, out) {
        s = tolower(h); sub(/^0x/, "", s); n = length(s); out = ""
        for (i = n; i >= 1; i--) {
            d = index("0123456789abcdef", substr(s, i, 1)) - 1
            for (b = 0; b < 4; b++) if (int(d / (2 ^ b)) % 2 == 1) out = out " " ((n - i) * 4 + b)
        }
        return out
    }
    function topo(c, campo,    f, v) {
        f = syscpu "/cpu" c "/topology/" campo; v = "?"
        if ((getline v < f) <= 0) v = "?"
        close(f); return v
    }
    {
        match($0, /task +[0-9]+/);  tk = substr($0, RSTART, RLENGTH); sub(/task +/, "", tk)
        match($0, /mask 0x[0-9a-fA-F]+/); mk = substr($0, RSTART + 5, RLENGTH - 5)
        if (tk in visto) next
        visto[tk] = 1; ntk++
        nc = split(hex2cpus(mk), cs, " ")
        split("", propio); nn = 0; txt = ""
        for (i = 1; i <= nc; i++) {
            c = cs[i]; pk = topo(c, "physical_package_id"); key = pk ":" topo(c, "core_id")
            txt = txt " " c
            if (!(key in propio)) { propio[key] = 1; nn++ }
            if (key in dueno) { if (dueno[key] != tk) compartido++ }
            else { dueno[key] = tk; porsock[pk]++; total++ }
        }
        printf "tarea %s mascara %s cpus_logicos%s nucleos_fisicos=%d\n", tk, mk, txt, nn
        if (nn != cpt) malcpt++
    }
    END {
        ok = (ntk == nt && compartido == 0 && malcpt == 0 && total == nt * cpt) ? 1 : 0
        s = ""
        for (k = 0; k <= 64; k++) if (k in porsock) s = s (s == "" ? "" : ",") "s" k ":" porsock[k]
        printf "afinidad_ok=%d\nprocesos_con_mascara=%d\nnucleos_fisicos_distintos=%d\nnucleos_compartidos=%d\nnucleos_por_socket=%s\n",
               ok, ntk, total + 0, compartido + 0, s
    }' > afinidad.txt
AF_OK=$(grep '^afinidad_ok=' afinidad.txt | cut -d= -f2)
SOCKS=$(grep '^nucleos_por_socket=' afinidad.txt | cut -d= -f2)

# ------------------------------------------------------------- 5b. donde corrieron los hilos
# Por muestra (misma hora), union de los nucleos fisicos de todos los hilos del job.
awk -v syscpu="$SYSCPU" '
    function topo(c, campo,    f, v) {
        f = syscpu "/cpu" c "/topology/" campo; v = "?"
        if ((getline v < f) <= 0) v = "?"
        close(f); return v
    }
    $1 == "muestra" {
        t = $2; mu[t] = 1; nh = 0
        for (i = 4; i <= NF; i++) {
            nh++; key = topo($i, "physical_package_id") ":" topo($i, "core_id")
            if (!((t SUBSEP key) in vis)) { vis[t, key] = 1; nuc[t]++ }
        }
        hil[t] += nh
    }
    END {
        ns = 0; mx = 0; mh = 0
        for (t in mu) { ns++; if (nuc[t] > mx) mx = nuc[t]; if (hil[t] > mh) mh = hil[t] }
        printf "hilos_muestras=%d\nhilos_totales=%d\nnucleos_fisicos_hilos=%d\n", ns, mh, mx
    }' hilos_crudo.txt > hilos.txt
HMUES=$(grep '^hilos_muestras=' hilos.txt | cut -d= -f2)
HNUC=$(grep '^nucleos_fisicos_hilos=' hilos.txt | cut -d= -f2)

# ------------------------------------------------------------- 6. carga ajena en los sockets del job
# En las muestras "durante" se restan los hilos propios (procesos x hilos).
AJENA=$(awk -v socks="$SOCKS" -v propios=$((NTASKS * NHILOS)) '
    BEGIN { n = split(socks, a, ","); for (i = 1; i <= n; i++) { split(a[i], b, ":"); js[b[1]] = 1 } }
    { oc = 0; tt = 0
      for (i = 3; i <= NF; i++) if ($i ~ /^s[0-9]+=/) {
          split($i, x, "="); split(x[2], y, "/")
          if (x[1] in js) { oc += y[1]; tt += y[2] } }
      if ($1 == "durante") oc -= propios
      if (oc < 0) oc = 0
      if (oc > mx) mx = oc
      suma += oc; nl++; T = tt }
    END { printf "carga_ajena_max=%d/%d\ncarga_ajena_media=%.1f\n", mx, T, (nl ? suma / nl : 0) }' carga.txt)

# ------------------------------------------------------------- 7. verificacion
FIN=$(grep -c "General timing and accounting" OUTCAR 2>/dev/null); FIN=${FIN:-0}
CONV=$(grep -c "aborting loop because EDIFF is reached" OUTCAR 2>/dev/null); CONV=${CONV:-0}
NSCF=$(grep -c "LOOP:" OUTCAR 2>/dev/null); NSCF=${NSCF:-0}
NMASK=$(grep -o "mask 0x[0-9a-fA-F]*" vasp.err 2>/dev/null | sort -u | wc -l)

{
    echo "fecha_fin=$(date '+%Y-%m-%d %H:%M:%S')"
    echo "wall_s=$WALL"
    echo "rc_srun=$RC"
    echo "ntasks=$NTASKS"
    echo "cpus_per_task=$CPT"
    echo "hilos=$NHILOS"
    echo "alloc_cpus=${ALLOC:-?}"
    echo "mascaras_distintas=$NMASK"
    grep -E '^(afinidad_ok|nucleos_fisicos_distintos|nucleos_compartidos|nucleos_por_socket)=' afinidad.txt
    cat hilos.txt
    echo "$AJENA"
    echo "fin_normal_vasp=$FIN"
    echo "convergido=$CONV"
    echo "n_scf=$NSCF"
} > tiempo.txt

echo "---- afinidad (srun --cpu-bind=verbose) ----"
cat afinidad.txt
echo "---- hilos ----"
cat hilos.txt
echo "---- carga del nodo ----"
cat carga.txt
echo "---- tiempo ----"
cat tiempo.txt

[ "$RC" -eq 0 ]            || { tail -20 vasp.out; tail -20 vasp.err; terminar 3 "FALLA_srun_rc_$RC"; }
[ "$FIN" -ge 1 ]           || { tail -20 vasp.out; terminar 3 "FALLA_VASP_no_termino_normalmente"; }
[ "$NSCF" -ge 1 ]          || terminar 3 "FALLA_sin_iteraciones_SCF"
[ "${AF_OK:-0}" = "1" ]    || terminar 4 "FALLA_afinidad_(ver_afinidad.txt)"
if [ "$NHILOS" -gt 1 ] && [ "${HMUES:-0}" -gt 0 ] && [ "${HNUC:-0}" -lt $((NTASKS * NHILOS)) ]; then
    terminar 4 "FALLA_hilos_en_${HNUC}_nucleos_fisicos_(se_esperaban_$((NTASKS * NHILOS)))"
fi
if [ -n "${NELM_DIAG:-}" ]; then
    terminar 0 "OK"
fi
[ "$CONV" -ge 1 ]          || terminar 5 "FALLA_SCF_no_convergio"
terminar 0 "OK"
