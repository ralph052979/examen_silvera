#!/bin/bash
# =============================================================================
#  empaquetar_datos.sh - junta en un solo archivo comprimido lo que hace falta para el
#  analisis final y el reporte:
#    * logs/ (registros de Slurm), datos/ (tablas y evidencia), reporte/figuras/
#    * sacct con los tiempos de cola de todos los jobs (para la D.1) y el historial de Git
#    * de cada corrida, solo los archivos pequenos: tiempo.txt, entorno.txt, afinidad.txt,
#      hilos.txt, progreso.txt, geometria_pasos.txt, carga.txt, INCAR, KPOINTS, OSZICAR y OUTCAR
#  No incluye el POTCAR (licencia de VASP) ni archivos grandes (WAVECAR, CHGCAR, vasprun.xml).
#  No usa VASP ni envia jobs (sacct solo consulta).
#
#  Uso (desde la carpeta del repositorio):   bash slurm/empaquetar_datos.sh
#  Deja el archivo en la carpeta de arriba del repositorio:  ../examen_datos.tar.gz
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
SALIDA="$(dirname "$REPO")/examen_datos.tar.gz"
EXTRA=logs/empaque
falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }

echo "=== 1/4  Git ==="
git add slurm/empaquetar_datos.sh analisis/parteA2_hilos.py README.md 2>/dev/null
if ! git diff --cached --quiet; then
    git commit -q -m "Script para empaquetar los datos del analisis final; aviso corregido en parteA2_hilos.py" \
        || falla "git commit fallo"
    echo "  [OK] commit $(git rev-parse --short HEAD)"
else
    echo "  [OK] nada nuevo que guardar"
fi

echo "=== 2/4  Consultas (sacct y Git) ==="
mkdir -p "$EXTRA"
sacct -u "$(id -un)" -S 2026-11-16 -X -P \
      --format=JobID,JobName%40,Partition,AllocCPUS,ReqMem,Timelimit,Submit,Eligible,Start,End,Elapsed,CPUTimeRAW,State,NodeList \
      > "$EXTRA/sacct_colas.txt" 2>&1 || falla "sacct fallo: $(head -2 "$EXTRA/sacct_colas.txt")"
sacct -u "$(id -un)" -S 2026-11-16 -P \
      --format=JobID,JobName%40,AllocCPUS,Elapsed,MaxRSS,AveRSS,MaxRSSTask,NTasks,State \
      > "$EXTRA/sacct_pasos.txt" 2>&1
git log --date=iso --format='%h %ad %s' > "$EXTRA/git_log.txt" 2>&1
git status --short > "$EXTRA/git_status.txt" 2>&1
echo "  [OK] $(($(wc -l < "$EXTRA/sacct_colas.txt") - 1)) jobs en sacct; $(wc -l < "$EXTRA/git_log.txt") commits en Git"

echo "=== 3/4  Archivos de las corridas ==="
LISTA="$EXTRA/lista.txt"
find corridas -mindepth 2 -maxdepth 2 -type f \( -name tiempo.txt -o -name entorno.txt -o -name afinidad.txt \
     -o -name hilos.txt -o -name progreso.txt -o -name geometria_pasos.txt -o -name carga.txt -o -name INCAR \
     -o -name KPOINTS -o -name OSZICAR -o -name OUTCAR \) | sort > "$LISTA"
N=$(cut -d/ -f2 "$LISTA" | sort -u | wc -l)
[ "$N" -gt 0 ] || falla "no encontre corridas en corridas/"
echo "  [OK] $N corridas, $(wc -l < "$LISTA") archivos"

echo "=== 4/4  Archivo comprimido ==="
rm -f "$SALIDA"
FIGS=""
[ -d reporte/figuras ] && FIGS=reporte/figuras
tar --exclude='POTCAR*' -czf "$SALIDA" logs datos $FIGS -T "$LISTA" || falla "tar fallo"
if tar -tzf "$SALIDA" | grep -q 'POTCAR'; then
    rm -f "$SALIDA"
    falla "el archivo incluia un POTCAR y se borro"
fi
echo "  [OK] $SALIDA  ($(du -h "$SALIDA" | cut -f1))"
[ "$(du -m "$SALIDA" | cut -f1)" -le 30 ] || echo "  [AVISO] pesa mas de 30 MB: si no se puede subir al chat, avisame"
echo
echo "Bajalo con WinSCP (esta en $(dirname "$REPO")) y subelo al chat."
