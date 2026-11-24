#!/bin/bash
# =============================================================================
#  precisar_mpirun.sh - guarda en Git la precision sobre el lanzador MPI en README.md y en el reporte:
#  mpirun SIN --map-by ni --bind-to ata todos los procesos al nucleo 0; la forma que indico el docente
#  (mpirun --map-by slot:PE=N --bind-to core) fija el mapeo de forma explicita. Con la direccion del
#  repositorio, ademas lo sube (slurm/publicar.sh).
#
#  Uso (desde la carpeta del repositorio):
#      bash slurm/precisar_mpirun.sh                                                  # solo el commit
#      bash slurm/precisar_mpirun.sh https://github.com/USUARIO/examen_silvera.git     # commit y subida
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }
ok()    { echo "  [OK] $1"; }
paso()  { echo; echo "=== $1 ==="; }
URL="${1:-}"
PROHIBIDOS='(^|/)(POTCAR|WAVECAR|CHGCAR|CHG|vasprun\.xml|PROCAR|DOSCAR|EIGENVAL)$'

paso "1/2  Revision de los cambios"
grep -q 'sin `--map-by` ni `--bind-to`' README.md || falla "README.md es una version anterior: descomprime de nuevo el paquete"
grep -q 'cod{mpirun} sin' reporte/reporte.tex || falla "reporte.tex es una version anterior: descomprime de nuevo el paquete"
[ -s reporte/reporte.pdf ] || falla "falta reporte/reporte.pdf"
ok "README.md y reporte.tex traen la frase corregida"

paso "2/2  Guardar en Git"
git add README.md reporte/reporte.tex reporte/reporte.pdf slurm/precisar_mpirun.sh || falla "git add fallo"
MAL=$(git diff --cached --name-only | grep -E "$PROHIBIDOS")
if [ -n "$MAL" ]; then
    git reset -q
    echo "$MAL"
    falla "se iban a agregar archivos prohibidos (arriba); no se agrego nada"
fi
if [ -n "$(git diff --cached --name-only)" ]; then
    git commit -q -m "Precisa el lanzador: mpirun sin --map-by ni --bind-to ata todo al nucleo 0 (README y reporte)" \
        || falla "git commit fallo"
    ok "commit $(git rev-parse --short HEAD)"
else
    ok "nada nuevo que guardar (ya estaba en Git)"
fi

if [ -z "$URL" ]; then
    echo
    echo "  Para subirlo a GitHub:"
    echo "      bash slurm/precisar_mpirun.sh https://github.com/TU_USUARIO/examen_silvera.git"
    exit 0
fi
exec bash slurm/publicar.sh "$URL"
