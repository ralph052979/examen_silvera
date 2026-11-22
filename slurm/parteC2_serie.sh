#!/bin/bash
#SBATCH --partition=phys_fast
#SBATCH --nodes=1
#SBATCH --output=logs/%x_%j.out
#SBATCH --signal=B:USR1@150
# =============================================================================
#  parteC2_serie.sh - Tarea C.2: una serie de ENCUT en un solo job de 1 nucleo.
#
#  Corre slurm/vasp_job.sh una vez por cada ENCUT de la variable ENCUTS (separados por
#  ':'), con CASO=encut/E<ENCUT> y la carpeta corridas/<job>_<id>_E<ENCUT>/, de modo que
#  cada corrida tiene las mismas verificaciones que las de las Partes A y B. Las corridas
#  van una tras otra, en el orden dado, en el mismo nucleo.
#  Si Slurm avisa que se acaba el tiempo (USR1), se detiene la corrida en curso (con sus
#  verificaciones) y no se empiezan las siguientes.
#
#  Lo envia slurm/parteC_enviar.sh con:  ENCUTS, GIT_HASH, MAX_PRIMERA_ITER, BLAS_HILOS.
# =============================================================================
set -u
REPO="${SLURM_SUBMIT_DIR:?este script se envia con sbatch desde la carpeta del repositorio}"
cd "$REPO" || exit 1
PARAR=0; HIJO=""
avisar() { PARAR=1; [ -n "$HIJO" ] && kill -USR1 "$HIJO" 2>/dev/null; }
trap avisar USR1
LISTA=$(echo "${ENCUTS:?falta la variable ENCUTS}" | tr ':' ' ')
echo "serie_encut: $LISTA"
for E in $LISTA; do
    if [ "$PARAR" = 1 ]; then
        echo "SERIE E$E: OMITIDA (se acabo el tiempo del job)"
        continue
    fi
    CASO="encut/E$E" SUFIJO="_E$E" HILOS=1 bash "$REPO/slurm/vasp_job.sh" &
    HIJO=$!
    while :; do
        wait "$HIJO"; RC=$?
        kill -0 "$HIJO" 2>/dev/null || break
    done
    HIJO=""
    echo "SERIE E$E: codigo $RC, $(grep '^estado=' "corridas/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_E$E/tiempo.txt" 2>/dev/null | tail -1)"
done
