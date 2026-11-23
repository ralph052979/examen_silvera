#!/bin/bash
# =============================================================================
#  publicar.sh - deja el repositorio listo para GitHub o GitLab y lo sube con todo su
#  historial de Git.
#   1. Agrega a Git, de cada corrida, solo los archivos pequenos (tiempo.txt, entorno.txt,
#      afinidad.txt, hilos.txt, progreso.txt, INCAR, KPOINTS, OSZICAR, OUTCAR, ...) y las
#      listas de jobs de logs/, para que el analisis se pueda repetir desde un clon.
#      Nunca agrega el POTCAR (licencia de VASP) ni archivos grandes (.gitignore).
#   2. Si se da la direccion de un repositorio VACIO de GitHub o GitLab, lo sube (git push).
#
#  Uso (desde la carpeta del repositorio):
#      bash slurm/publicar.sh                                              # solo el paso 1
#      bash slurm/publicar.sh https://github.com/USUARIO/examen_silvera.git # paso 1 y subida
#  Al subir, Git pide el usuario de GitHub y una "contrasena": ahi se pega el token
#  (no se ve mientras se escribe).
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }
ok()    { echo "  [OK] $1"; }
paso()  { echo; echo "=== $1 ==="; }
URL="${1:-}"
PROHIBIDOS='(^|/)(POTCAR|WAVECAR|CHGCAR|CHG|vasprun\.xml|PROCAR|DOSCAR|EIGENVAL)$'

paso "1/3  Revision"
grep -q '^!corridas/\*/OUTCAR' .gitignore || falla ".gitignore es una version anterior: descomprime de nuevo el paquete"
[ -z "$(git ls-files | grep -E "$PROHIBIDOS")" ] || { git ls-files | grep -E "$PROHIBIDOS"; falla "hay archivos prohibidos en Git (arriba)"; }
ok "ningun POTCAR ni archivo grande de VASP en Git"
[ -n "$(git config user.name)" ] || falla "Git no tiene tu nombre: git config --global user.name \"Rafael Silvera Altamiranda\""

paso "2/3  Archivos pequenos de las corridas y listas de jobs"
git add .gitignore README.md slurm/publicar.sh corridas logs 2>/dev/null || falla "git add fallo"
MAL=$(git diff --cached --name-only | grep -E "$PROHIBIDOS")
if [ -n "$MAL" ]; then
    git reset -q
    echo "$MAL"
    falla "se iban a agregar archivos prohibidos (arriba); no se agrego nada"
fi
N=$(git diff --cached --name-only | wc -l)
if [ "$N" -gt 0 ]; then
    TAM=$(git diff --cached --name-only -z | xargs -0 du -cb 2>/dev/null | tail -1 | awk '{printf "%.1f MB", $1 / 1048576}')
    git commit -q -m "Archivos pequenos de cada corrida y listas de jobs (analisis reproducible desde un clon)" \
        || falla "git commit fallo"
    ok "commit $(git rev-parse --short HEAD): $N archivos ($TAM)"
else
    ok "nada nuevo que guardar (ya estaba en Git)"
fi
ok "corridas en Git: $(git ls-files corridas | cut -d/ -f2 | sort -u | wc -l); tamano del historial: $(du -sh .git | cut -f1)"
[ -z "$(git status --porcelain -- analisis teoria reporte datos slurm README.md entradas)" ] \
    || falla "hay cambios sin guardar en Git (git status)"

paso "3/3  Subir a GitHub o GitLab"
if [ -z "$URL" ]; then
    echo "  Falta la direccion del repositorio. Crea en GitHub un repositorio VACIO llamado examen_silvera"
    echo "  (sin README, sin .gitignore y sin licencia) y corre:"
    echo "      bash slurm/publicar.sh https://github.com/TU_USUARIO/examen_silvera.git"
    exit 0
fi
case "$URL" in
    *github.com*|*gitlab*) ;;
    *) echo "  [AVISO] la direccion no es de GitHub ni de GitLab: $URL" ;;
esac
if git remote | grep -qx origin; then
    git remote set-url origin "$URL" || falla "no se pudo cambiar la direccion de origin"
else
    git remote add origin "$URL" || falla "no se pudo agregar la direccion"
fi
echo "  Git va a pedir el usuario de GitHub y la contrasena: en la contrasena pega el token"
echo "  (no se ve mientras se escribe) y presiona Enter."
git push -u origin HEAD:main || falla "git push fallo: revisa el usuario, el token y que el repositorio de GitHub este vacio"
WEB=$(echo "$URL" | sed -e 's/\.git$//')
echo
echo "================================================================"
echo " REPOSITORIO SUBIDO: $WEB"
echo " Ese es el enlace que se entrega junto con reporte/reporte.pdf"
echo "================================================================"
