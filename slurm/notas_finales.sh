#!/bin/bash
# =============================================================================
#  notas_finales.sh - guarda en Git la nota que separa, en el registro de sacct del enunciado
#  (datos/sacct_registro_examen.txt, desde el 2026-09-07), los jobs del examen (JobID >= 528789)
#  de los jobs de investigacion anteriores, las versiones de Python en el README y el reporte con
#  la seccion de discusion. Con la direccion del repositorio, ademas lo sube (slurm/publicar.sh).
#
#  Uso (desde la carpeta del repositorio):
#      bash slurm/notas_finales.sh                                                  # solo el commit
#      bash slurm/notas_finales.sh https://github.com/USUARIO/examen_silvera.git     # commit y subida
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }
ok()    { echo "  [OK] $1"; }
paso()  { echo; echo "=== $1 ==="; }
URL="${1:-}"
REG=datos/sacct_registro_examen.txt
PRIMER_JOB=528789
ESPERADO="370 11682.7 166 148.5"
PROHIBIDOS='(^|/)(POTCAR|WAVECAR|CHGCAR|CHG|vasprun\.xml|PROCAR|DOSCAR|EIGENVAL)$'

paso "1/2  Revision de la nota"
[ -s "$REG" ] || falla "no existe $REG"
[ -s datos/sacct_registro_LEEME.md ] || falla "falta datos/sacct_registro_LEEME.md: descomprime de nuevo el paquete"
grep -q 'sacct_registro_LEEME' README.md || falla "README.md es una version anterior: descomprime de nuevo el paquete"
grep -q 'section{Discusi' reporte/reporte.tex || falla "reporte.tex es una version anterior: descomprime de nuevo el paquete"
[ -s reporte/reporte.pdf ] || falla "falta reporte/reporte.pdf"
MEDIDO=$(awk -v p="$PRIMER_JOB" '$1~/^[0-9]/ {h=$6/3600; if ($1+0<p) {n1++; h1+=h} else {n2++; h2+=h}}
    END {printf "%d %.1f %d %.1f", n1, h1, n2, h2}' "$REG")
set -- $MEDIDO
echo "  registro del enunciado: $(grep -c '^ *[0-9]' "$REG") jobs"
echo "    antes del job $PRIMER_JOB (investigacion): $1 jobs, $2 h-nucleo"
echo "    desde el job $PRIMER_JOB (examen):         $3 jobs, $4 h-nucleo (de 250)"
[ "$MEDIDO" = "$ESPERADO" ] || falla "los numeros no coinciden con la nota ($ESPERADO)"
ok "la nota coincide con el registro de sacct"

paso "2/2  Guardar en Git"
git add README.md datos/sacct_registro_LEEME.md reporte/reporte.tex reporte/reporte.pdf slurm/notas_finales.sh \
    || falla "git add fallo"
MAL=$(git diff --cached --name-only | grep -E "$PROHIBIDOS")
if [ -n "$MAL" ]; then
    git reset -q
    echo "$MAL"
    falla "se iban a agregar archivos prohibidos (arriba); no se agrego nada"
fi
if [ -n "$(git diff --cached --name-only)" ]; then
    git commit -q -m "Nota del registro de sacct (jobs del examen desde 528789), versiones de Python y seccion de discusion" \
        || falla "git commit fallo"
    ok "commit $(git rev-parse --short HEAD)"
else
    ok "nada nuevo que guardar (ya estaba en Git)"
fi

if [ -z "$URL" ]; then
    echo
    echo "  Para subirlo a GitHub:"
    echo "      bash slurm/notas_finales.sh https://github.com/TU_USUARIO/examen_silvera.git"
    exit 0
fi
exec bash slurm/publicar.sh "$URL"
