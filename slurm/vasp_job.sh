#!/bin/bash
#SBATCH --partition=phys_fast
#SBATCH --nodes=1
#SBATCH --output=logs/%x_%j.out
#SBATCH --signal=B:USR1@90
#
# =============================================================================
#  vasp_job.sh - ejecuta UN calculo VASP en su propia carpeta
#                corridas/<nombre_del_job>_<jobid>/
#
#  No se envia a mano: lo envian slurm/enviar.sh (Tarea 0),
#  slurm/parteA1_enviar.sh (A.1) y los scripts parteA2_*.sh (A.2). Ellos fijan
#  --job-name, --ntasks, --cpus-per-task, --hint o --threads-per-core, --mem,
#  --time y las variables de abajo.
#
#  Variables de entrada (sbatch --export=ALL,...):
#    CASO         K1..K7     -> entradas/candidatos/<CASO>/  (INCAR KPOINTS POSCAR)
#                 caso_base  -> entradas/caso_base/
#    GIT_HASH     (opcional) commit del repositorio al momento de enviar
#    PASO_NTASKS  (opcional) procesos del paso srun   (por defecto, los del job)
#    PASO_CPT     (opcional) CPUs por proceso del paso (por defecto, las del job)
#    PASO_TPC     (opcional) srun --threads-per-core=PASO_TPC
#    HILOS        (opcional) hilos OpenMP por proceso  (por defecto PASO_CPT)
#    BLAS_HILOS   (opcional) OPENBLAS_NUM_THREADS (1 = OpenBLAS sin hilos propios)
#    ESPERA_OMP   (opcional) OMP_WAIT_POLICY (PASSIVE o ACTIVE; vacio = la de libgomp)
#    NELM_DIAG    (opcional) fija NELM: corrida corta de diagnostico (no exige convergencia)
#    MAX_PRIMERA_ITER (opcional) segundos: si VASP no termino ninguna iteracion SCF
#                 en ese tiempo, se detiene la corrida (ahorra horas-nucleo)
#    REQUIERE_OK / REQUIERE_FALLA (opcional) numeros de job separados por ':'.
#                 La corrida se omite (en segundos) si alguno de REQUIERE_OK no
#                 termino OK, o si alguno de REQUIERE_FALLA si termino OK.
#    PILA_HILOS   (opcional) n = en la muestra n de hilos, foto de las pilas de todos
#                 los hilos con gdb (pila_hilos.txt; solo diagnostico)
#
#  Configuracion de Granado verificada (ver README.md):
#    * --hint=nomultithread    : 1 proceso por nucleo FISICO; AllocCPUS = ntasks
#    * srun --mpi=pmix         : unico lanzador que respeta la asignacion de Slurm
#    * lib/ en LD_LIBRARY_PATH : libfftw3_omp.so.3 no existe en los nodos de computo
#    * tiempo con date         : /usr/bin/time no existe en los nodos
#    * NO fijar PMIX_MCA_gds   : rompe MPI_Init. Los avisos "PMIX ERROR: BAD-PARAM"
#                                son inofensivos.
#
#  Verificaciones de cada corrida (quedan en tiempo.txt, afinidad.txt e hilos.txt):
#    * afinidad: cada proceso tiene exactamente HILOS nucleos FISICOS en su mascara
#      y ningun nucleo fisico se comparte entre procesos (se mira al arrancar:
#      si falla, la corrida se detiene en segundos);
#    * hilos (si HILOS > 1): cada 10 s se lee /proc/<pid>/task/*/stat de los
#      procesos de VASP (el pid lo escribe srun en vasp.err). Se mide en que
#      nucleos fisicos estan los hilos y cuantos hilos estan ACTIVOS (>= 50 % de
#      un CPU). Si hay mas hilos activos que nucleos (p. ej. un hilo propio de
#      OpenBLAS), la corrida se marca como sobresuscrita;
#    * VASP termino normalmente y el SCF convergio (salvo con NELM_DIAG);
#    * carga ajena: CPUs ocupados por otros usuarios en los sockets del job.
#  Si Slurm avisa que se acaba el tiempo (USR1, 90 s antes), se detiene VASP y se
#  completan igual las verificaciones.
# =============================================================================

set -u
REPO="${SLURM_SUBMIT_DIR:?este script se envia con sbatch desde la carpeta del repositorio}"
cd "$REPO" || exit 1
SYSCPU="${SYSCPU:-/sys/devices/system/cpu}"    # solo se cambia en pruebas
PROCSTAT="${PROCSTAT:-/proc/stat}"              # solo se cambia en pruebas
PROCPID="${PROCPID:-/proc}"                     # solo se cambia en pruebas
MON_INTERVALO="${MON_INTERVALO:-10}"            # s entre muestras de hilos
MON_INICIO="${MON_INICIO:-5}"                   # s antes de la primera muestra

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

# Corridas condicionales: se decide con el estado de jobs anteriores de la misma cadena.
estado_de() {
    local d
    d=$(ls -d "$REPO"/corridas/*_"$1" 2>/dev/null | head -1)
    [ -n "$d" ] && grep '^estado=' "$d/tiempo.txt" 2>/dev/null | tail -1 | cut -d= -f2
}
omitir() {
    echo "estado=OMITIDA_$1" > "$RUN/tiempo.txt"
    echo "RESULTADO: OMITIDA_$1"
    exit 0
}
for id in $(echo "${REQUIERE_OK:-}" | tr ':' ' '); do
    e=$(estado_de "$id")
    [ "$e" = "OK" ] || omitir "el_job_${id}_no_termino_OK_(${e:-sin_estado})"
done
for id in $(echo "${REQUIERE_FALLA:-}" | tr ':' ' '); do
    e=$(estado_de "$id")
    [ "$e" = "OK" ] && omitir "el_job_${id}_ya_termino_OK"
done

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
[ -n "${PASO_TPC:-}" ]    && GEOM="$GEOM --threads-per-core=$PASO_TPC"
# Solo este script decide los hilos: se borra lo que venga del entorno del login.
unset OPENBLAS_NUM_THREADS GOTO_NUM_THREADS OMP_WAIT_POLICY GOMP_SPINCOUNT OMP_DYNAMIC \
      OMP_NESTED OMP_MAX_ACTIVE_LEVELS MKL_NUM_THREADS
export OMP_NUM_THREADS="$NHILOS"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
[ -n "${BLAS_HILOS:-}" ] && export OPENBLAS_NUM_THREADS="$BLAS_HILOS"
[ -n "${ESPERA_OMP:-}" ] && export OMP_WAIT_POLICY="$ESPERA_OMP"
# VASP compilado con gfortran guarda OSZICAR y OUTCAR en un bufer hasta el final; la
# salida estandar (vasp.out) sin bufer muestra cada iteracion SCF en el momento (DAV: n).
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
ulimit -s unlimited 2>/dev/null || echo "aviso: no se pudo fijar ulimit -s unlimited"

cd "$RUN" || exit 1
VASP_BIN="$(command -v vasp_std)" || terminar 2 "FALLA_vasp_std_no_encontrado"
CPU_MODELO=$(grep -m1 "model name" /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ *//')
ldd "$VASP_BIN" > ldd.txt 2>&1
# OpenBLAS: con hilos de OpenMP (enlaza libgomp) o con hilos propios (pthreads)?
LIBOB=$(awk '/libopenblas/ && $2 == "=>" {print $3; exit}' ldd.txt)
OB_OMP="?"; OB_VER="?"
if [ -n "$LIBOB" ] && [ -f "$LIBOB" ]; then
    if ldd "$LIBOB" 2>/dev/null | grep -q 'libgomp'; then OB_OMP="si"; else OB_OMP="no"; fi
    OB_VER=$(LC_ALL=C grep -a -o -E -m1 'OpenBLAS [0-9]+\.[0-9]+\.[0-9]+[A-Za-z_ ]{0,60}' "$LIBOB" 2>/dev/null | head -1 | tr -s ' ' | sed 's/ *$//')
fi
FFTLIB=$(LC_ALL=C grep -a -c -i 'fftlib' "$VASP_BIN" 2>/dev/null)

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
    echo "paso_threads_per_core=${PASO_TPC:-no}"
    echo "omp_num_threads=$OMP_NUM_THREADS"
    echo "openblas_num_threads=${OPENBLAS_NUM_THREADS:-no_fijado}"
    echo "omp_wait_policy=${OMP_WAIT_POLICY:-no_fijado}"
    echo "openblas_lib=${LIBOB:-?}"
    echo "openblas_con_openmp=$OB_OMP"
    echo "openblas_version=${OB_VER:-?}"
    echo "cadena_fftlib_en_vasp=${FFTLIB:-?}"
    echo "max_primera_iter_s=${MAX_PRIMERA_ITER:-no}"
    echo "nelm_diag=${NELM_DIAG:-no}"
    echo "git_hash=${GIT_HASH:-desconocido}"
    echo "vasp_bin=$VASP_BIN"
    echo "fecha_inicio=$(date '+%Y-%m-%d %H:%M:%S')"
} > entorno.txt
module list > modulos.txt 2>&1
env | grep -E '^(SLURM_|OMP_|GOMP_|OPENBLAS_|GOTO_|MKL_|LD_LIBRARY_PATH)' | sort > variables.txt
scontrol show job -d "$SLURM_JOB_ID" > slurm_job.txt 2>&1
ALLOC=$(grep -o "NumCPUs=[0-9]*" slurm_job.txt | head -1 | cut -d= -f2)
md5sum INCAR KPOINTS POSCAR POTCAR > entradas.md5
for d in "$SYSCPU"/cpu[0-9]*; do
    echo "${d##*cpu} $(cat "$d/topology/physical_package_id" 2>/dev/null || echo '?') $(cat "$d/topology/core_id" 2>/dev/null || echo '?')"
done | sort -n > topologia.txt     # cpu socket nucleo

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

# ------------------------------------------------------------- 4. funciones de verificacion
# Afinidad: srun --cpu-bind=verbose escribe en vasp.err una linea por proceso:
#   cpu-bind=MASK - <nodo>, task  3  3 [12345]: mask 0x8000000008 set
afinidad() {
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
        }'
}
pids_vasp() {   # pid de cada proceso de VASP, segun las lineas de srun en vasp.err
    grep -oE "task +[0-9]+ +[0-9]+ \[[0-9]+\]: mask" vasp.err 2>/dev/null | grep -oE "\[[0-9]+\]" | tr -d '[]' | sort -un
}
detener() {   # detener <motivo>: deja el motivo y para VASP
    [ -s .motivo ] || echo "$1" > .motivo
    kill -TERM "$SRUN_PID" 2>/dev/null
    local i=0
    while kill -0 "$SRUN_PID" 2>/dev/null && [ "$i" -lt 30 ]; do sleep 1; i=$((i + 1)); done
    if kill -0 "$SRUN_PID" 2>/dev/null; then
        for q in $(pids_vasp); do kill -KILL "$q" 2>/dev/null; done
        kill -KILL "$SRUN_PID" 2>/dev/null
    fi
}

iteraciones() {   # iteraciones SCF que VASP ya escribio en la salida estandar
    local n
    n=$(grep -cE '^(DAV|RMM|DIA|SDA|CG) *: +[0-9]' vasp.out 2>/dev/null)
    echo "${n:-0}"
}

# Vigilante: (1) mascara correcta al arrancar; (2) anota en progreso.txt el segundo en que
# termina cada iteracion SCF; (3) si MAX_PRIMERA_ITER, detiene la corrida sin iteraciones.
vigilante() {
    local revisada=0 n vistas=0 k
    while kill -0 "$SRUN_PID" 2>/dev/null; do
        sleep 2
        if [ "$revisada" = 0 ]; then
            n=$(grep -cE "task +[0-9]+ +[0-9]+ \[[0-9]+\]: mask" vasp.err 2>/dev/null)
            if [ "${n:-0}" -ge "$NTASKS" ]; then
                revisada=1
                afinidad > .afinidad_inicio.txt
                if ! grep -q '^afinidad_ok=1' .afinidad_inicio.txt; then
                    detener afinidad
                    return
                fi
            fi
        fi
        n=$(iteraciones)
        if [ "$n" -gt "$vistas" ]; then
            for k in $(seq $((vistas + 1)) "$n"); do
                echo "iteracion $k $((SECONDS - INICIO_SRUN))" >> progreso.txt
            done
            vistas=$n
        fi
        if [ -n "${MAX_PRIMERA_ITER:-}" ] && [ "$n" -eq 0 ] \
           && [ $((SECONDS - INICIO_SRUN)) -ge "$MAX_PRIMERA_ITER" ]; then
            detener lento
            return
        fi
    done
}

# Monitor de hilos: cada MON_INTERVALO s, una ronda con "proc <pid> <uptime> tid:cpu:estado:ticks ...".
monitor_hilos() {
    local pids="" i=0 ronda=0 vivos pid lin
    while [ -z "$pids" ] && [ "$i" -lt 120 ]; do
        sleep 1; i=$((i + 1))
        pids=$(pids_vasp)
        kill -0 "$SRUN_PID" 2>/dev/null || break
    done
    [ -n "$pids" ] || { echo "sin_pids"; return; }
    echo "pids $(echo $pids)"
    sleep "$MON_INICIO"
    while :; do
        ronda=$((ronda + 1)); vivos=0
        echo "ronda $ronda"
        for pid in $pids; do
            [ -d "$PROCPID/$pid/task" ] || continue
            vivos=1
            lin=$(cat "$PROCPID/$pid"/task/*/stat 2>/dev/null | awk '{
                    s = $0; sub(/^[^)]*\) /, "", s); n = split(s, f, " ")
                    if (n >= 37) printf " %s:%s:%s:%d", $1, f[37], f[1], f[12] + f[13] }')
            echo "proc $pid $(cut -d' ' -f1 /proc/uptime)$lin"
        done
        [ "$vivos" = 1 ] || break
        if [ "$ronda" -eq "${PILA_HILOS:-0}" ]; then
            if command -v gdb >/dev/null 2>&1; then
                for pid in $pids; do
                    timeout 120 gdb -batch -nx -p "$pid" -ex "thread apply all bt 8" > pila_hilos.txt 2>&1 \
                        || echo "gdb termino con codigo $?" >> pila_hilos.txt
                    break
                done
            else
                echo "gdb no esta instalado en el nodo" > pila_hilos.txt
            fi
        fi
        sleep "$MON_INTERVALO"
    done
}

# ------------------------------------------------------------- 5. VASP
: > .motivo
: > hilos_crudo.txt
: > progreso.txt
INICIO_SRUN=$SECONDS
T0=$(date +%s.%N)
srun --mpi=pmix $GEOM --cpu-bind=verbose,cores vasp_std < /dev/null > vasp.out 2> vasp.err &
SRUN_PID=$!
trap 'detener limite_de_tiempo' USR1
vigilante &
VIG=$!
MONH=""
if [ "$NHILOS" -gt 1 ]; then
    monitor_hilos >> hilos_crudo.txt 2>/dev/null &
    MONH=$!
fi
RC=0
while :; do
    wait "$SRUN_PID"
    RC=$?
    # wait vuelve antes (codigo > 128) si llega USR1; entonces se espera de nuevo a srun
    if [ "$RC" -gt 128 ] && kill -0 "$SRUN_PID" 2>/dev/null; then continue; fi
    break
done
T1=$(date +%s.%N)
trap - USR1
for m in $MONITOR $MONH $VIG; do
    pkill -P "$m" 2>/dev/null
    kill "$m" 2>/dev/null
    wait "$m" 2>/dev/null
done
n=$(iteraciones); k=$(wc -l < progreso.txt)
while [ "$k" -lt "$n" ]; do k=$((k + 1)); echo "iteracion $k $((SECONDS - INICIO_SRUN))" >> progreso.txt; done
MOTIVO=$(cat .motivo 2>/dev/null)
carga fin >> carga.txt
WALL=$(awk -v a="$T0" -v b="$T1" 'BEGIN{printf "%.3f", b-a}')

# ------------------------------------------------------------- 6. afinidad por nucleo fisico
afinidad > afinidad.txt
AF_OK=$(grep '^afinidad_ok=' afinidad.txt | cut -d= -f2)
SOCKS=$(grep '^nucleos_por_socket=' afinidad.txt | cut -d= -f2)

# ------------------------------------------------------------- 6b. hilos: donde estan y cuantos trabajan
# Por ronda: nucleos fisicos de TODOS los hilos (donde quedaron atados) y, entre dos
# rondas, hilos ACTIVOS (>= 50 % de un CPU) y uso total de CPU de los procesos de VASP.
awk -v syscpu="$SYSCPU" -v hz="$(getconf CLK_TCK 2>/dev/null || echo 100)" -v lim=$((NTASKS * NHILOS)) '
    function topo(c, campo,    f, v) {
        f = syscpu "/cpu" c "/topology/" campo; v = "?"
        if ((getline v < f) <= 0) v = "?"
        close(f); return v
    }
    function cerrar(   k) {
        if (!abierta) return
        if (nh > mh) mh = nh
        if (nuc > mx) mx = nuc
        if (completa && npr > 0) { ni++; A[ni] = act; U[ni] = uso }
        abierta = 0
    }
    function mediana(v, n,    i, j, x, w) {
        if (n == 0) return 0
        for (i = 1; i <= n; i++) w[i] = v[i]
        for (i = 2; i <= n; i++) { x = w[i]; j = i - 1; while (j >= 1 && w[j] > x) { w[j + 1] = w[j]; j-- }; w[j + 1] = x }
        return (n % 2) ? w[(n + 1) / 2] : (w[n / 2] + w[n / 2 + 1]) / 2
    }
    BEGIN { if (hz + 0 <= 0) hz = 100 }
    $1 == "ronda" { cerrar(); abierta = 1; ns++; nh = 0; nuc = 0; act = 0; uso = 0; npr = 0; completa = 1; split("", vis); next }
    $1 == "proc" && abierta {
        pid = $2; t = $3 + 0; npr++
        valido = (pid in tprev) && (t > tprev[pid]); dt = t - tprev[pid]
        if (!valido) completa = 0
        for (i = 4; i <= NF; i++) {
            split($i, a, ":"); nh++
            key = topo(a[2], "physical_package_id") ":" topo(a[2], "core_id")
            if (!(key in vis)) { vis[key] = 1; nuc++ }
            if (valido && ((pid SUBSEP a[1]) in kprev)) {
                sh = (a[4] - kprev[pid, a[1]]) / (dt * hz)
                uso += sh; if (sh >= 0.5) act++
            }
            kprev[pid, a[1]] = a[4]
        }
        tprev[pid] = t
    }
    END {
        cerrar()
        am = mediana(A, ni); um = mediana(U, ni); amax = 0; umax = 0
        for (i = 1; i <= ni; i++) { if (A[i] > amax) amax = A[i]; if (U[i] > umax) umax = U[i] }
        sob = (ni >= 2 && (am > lim || um > lim + 0.5)) ? 1 : 0
        printf "hilos_muestras=%d\nhilos_totales=%d\nnucleos_fisicos_hilos=%d\n", ns, mh, mx
        printf "hilos_intervalos=%d\nhilos_activos_mediana=%g\nhilos_activos_max=%d\n", ni, am, amax
        printf "uso_cpu_mediana=%.2f\nuso_cpu_max=%.2f\nsobresuscripcion=%d\n", um, umax, sob
    }' hilos_crudo.txt > hilos.txt
HMUES=$(grep '^hilos_muestras=' hilos.txt | cut -d= -f2)
HNUC=$(grep '^nucleos_fisicos_hilos=' hilos.txt | cut -d= -f2)
HACT=$(grep '^hilos_activos_mediana=' hilos.txt | cut -d= -f2)
HUSO=$(grep '^uso_cpu_mediana=' hilos.txt | cut -d= -f2)
SOBRE=$(grep '^sobresuscripcion=' hilos.txt | cut -d= -f2)

# ------------------------------------------------------------- 7. carga ajena en los sockets del job
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

# ------------------------------------------------------------- 8. verificacion
FIN=$(grep -c "General timing and accounting" OUTCAR 2>/dev/null); FIN=${FIN:-0}
CONV=$(grep -c "aborting loop because EDIFF is reached" OUTCAR 2>/dev/null); CONV=${CONV:-0}
NSCF=$(grep -c "LOOP:" OUTCAR 2>/dev/null); NSCF=${NSCF:-0}
NOUT=$(iteraciones)
[ "$NOUT" -gt "$NSCF" ] && NSCF=$NOUT     # si VASP se detuvo, el OUTCAR puede quedar sin vaciar
NMASK=$(grep -o "mask 0x[0-9a-fA-F]*" vasp.err 2>/dev/null | sort -u | wc -l)

{
    echo "fecha_fin=$(date '+%Y-%m-%d %H:%M:%S')"
    echo "wall_s=$WALL"
    echo "rc_srun=$RC"
    echo "detenida_por=${MOTIVO:-no}"
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

LIM=$((NTASKS * NHILOS))
case "$MOTIVO" in
    afinidad)          cat afinidad.txt; terminar 4 "FALLA_afinidad_(detenida_al_arrancar,_ver_afinidad.txt)" ;;
    lento)             terminar 6 "FALLA_lento_(sin_iteraciones_SCF_en_${MAX_PRIMERA_ITER}_s)" ;;
    limite_de_tiempo)  terminar 6 "FALLA_limite_de_tiempo_(${NSCF}_iteraciones_SCF)" ;;
esac
[ "$RC" -eq 0 ]            || { tail -20 vasp.out; tail -20 vasp.err; terminar 3 "FALLA_srun_rc_$RC"; }
[ "$FIN" -ge 1 ]           || { tail -20 vasp.out; terminar 3 "FALLA_VASP_no_termino_normalmente"; }
[ "$NSCF" -ge 1 ]          || terminar 3 "FALLA_sin_iteraciones_SCF"
[ "${AF_OK:-0}" = "1" ]    || terminar 4 "FALLA_afinidad_(ver_afinidad.txt)"
if [ "$NHILOS" -gt 1 ] && [ "${HMUES:-0}" -gt 0 ] && [ "${HNUC:-0}" -lt "$LIM" ]; then
    terminar 4 "FALLA_hilos_en_${HNUC}_nucleos_fisicos_(se_esperaban_${LIM})"
fi
if [ "$NHILOS" -gt 1 ] && [ "${SOBRE:-0}" = "1" ]; then
    terminar 4 "FALLA_sobresuscripcion_(${HACT}_hilos_activos_y_uso_de_CPU_${HUSO}_para_${LIM}_nucleos)"
fi
if [ -n "${NELM_DIAG:-}" ]; then
    terminar 0 "OK"
fi
[ "$CONV" -ge 1 ]          || terminar 5 "FALLA_SCF_no_convergio"
terminar 0 "OK"
