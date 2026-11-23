#!/bin/bash
# =============================================================================
#  finalizar.sh - analisis final del examen. Corre, con los datos de corridas/ y datos/,
#  los scripts de analisis/ que generan las tablas y figuras de las Partes A.3, B, C y D,
#  las tablas del reporte (reporte/tablas/) y los calculos de la Parte I, y guarda todo en
#  Git junto con el reporte (reporte/reporte.tex y reporte.pdf).
#  No usa VASP ni envia jobs: solo Python en el login (sacct solo consulta).
#
#  Uso (desde la carpeta del repositorio):   bash slurm/finalizar.sh
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }
ok()    { echo "  [OK] $1"; }
paso()  { echo; echo "=== $1 ==="; }

paso "1/4  Archivos y Python"
for f in analisis/parteB_resumen.py analisis/parteA3_graficas.py analisis/parteB_graficas.py analisis/parteC_analisis.py \
         analisis/parteD_colas.py analisis/reporte_tablas.py teoria/parte1.md teoria/parte1_calculos.py \
         reporte/reporte.tex reporte/reporte.pdf datos/sacct_colas.txt README.md; do
    [ -s "$f" ] || falla "falta $f: el paquete no se descomprimio dentro de examen_silvera (ver la guia)"
done
CR=$(grep -l $'\r' slurm/*.sh analisis/*.py teoria/*.py teoria/*.md reporte/*.tex README.md 2>/dev/null)
[ -z "$CR" ] || falla "estos archivos tienen finales de linea de Windows: $CR (no los abras en Windows)"
for f in analisis/*.py teoria/*.py; do python3 -m py_compile "$f" 2>/dev/null || falla "error de sintaxis en $f"; done
rm -rf analisis/__pycache__ teoria/__pycache__
V=$(python3 -c "import matplotlib, numpy; print('matplotlib %s, numpy %s' % (matplotlib.__version__, numpy.__version__))" 2>/dev/null) \
    || falla "falta matplotlib o numpy: instala con  python3 -m pip install --user matplotlib numpy  y repite"
ok "$V"
[ -s logs/jobs_B1.txt ] && [ -s logs/jobs_C2.txt ] && [ -s logs/jobs_A1.txt ] || falla "faltan registros en logs/ (jobs_A1, jobs_B1, jobs_C2)"
ok "scripts, datos y registros en orden"

paso "2/4  Analisis (tablas y figuras)"
mkdir -p logs/final
correr() {   # correr <nombre> <script> : corre el script y guarda su salida en logs/final/<nombre>.txt
    python3 "$2" > "logs/final/$1.txt" 2>&1 || { tail -15 "logs/final/$1.txt"; falla "$2 termino con error"; }
}
correr B_resumen analisis/parteB_resumen.py
grep -q "VERIFICACION: las 40 corridas" logs/final/B_resumen.txt || { tail -8 logs/final/B_resumen.txt; falla "la verificacion de la Parte B no dio 40 corridas"; }
ok "Parte B: $(grep 'Combinacion optima' logs/final/B_resumen.txt | sed 's/ *(menor t medio)//')"
correr A3 analisis/parteA3_graficas.py
ok "A.3: $(grep -m1 '^MPI corrida' logs/final/A3.txt | awk '{print "Amdahl, f(MPI) = " $4 " +/- " $6}')"
correr B_graficas analisis/parteB_graficas.py
ok "B: $(grep -E 'pasa el limite|no cruza' logs/final/B_graficas.txt | sed 's/^ *//')"
correr C analisis/parteC_analisis.py
ok "C: $(grep -m1 '^24 ' logs/final/C.txt | awk '{print "p = 24: E(24) - E(1) = " $6 " eV = " $8 " eps_mach |E|"}')"
correr D analisis/parteD_colas.py
ok "D: $(grep -m1 '^Consumo' logs/final/D.txt)"
[ -s datos/sacct_registro_examen.txt ] && ok "registro de consumo del enunciado: datos/sacct_registro_examen.txt"
correr tablas analisis/reporte_tablas.py
ok "tablas del reporte: $(ls reporte/tablas/*.tex | wc -l) archivos en reporte/tablas/"
python3 teoria/parte1_calculos.py > teoria/parte1_salida.txt 2>&1 || falla "teoria/parte1_calculos.py termino con error"
ok "Parte I: teoria/parte1_salida.txt ($(grep -m1 'S(24)' teoria/parte1_salida.txt | sed 's/^ *//'))"
N=$(ls reporte/figuras/*.pdf 2>/dev/null | wc -l)
[ "$N" -ge 5 ] || falla "se esperaban al menos 5 figuras en reporte/figuras/ y hay $N"
ok "$N figuras en reporte/figuras/ (pdf y png)"

paso "3/4  Reporte"
PAG=""
command -v pdfinfo >/dev/null 2>&1 && PAG=$(pdfinfo reporte/reporte.pdf 2>/dev/null | awk '/^Pages:/ {print $2}')
ok "reporte/reporte.pdf ($(du -h reporte/reporte.pdf | cut -f1)${PAG:+, $PAG paginas}); fuente: reporte/reporte.tex"

paso "4/4  Git"
git add analisis teoria reporte datos slurm README.md .gitignore 2>/dev/null || falla "git add fallo"
if git diff --cached --quiet; then
    ok "nada nuevo que guardar (ya estaba en Git)"
else
    git commit -q -m "Analisis final (A.3, B, C y D), Parte I y reporte" || falla "git commit fallo"
    ok "commit $(git rev-parse --short HEAD)"
fi
[ -z "$(git status --porcelain -- analisis teoria reporte datos slurm README.md)" ] || {
    git status --short -- analisis teoria reporte datos slurm README.md | head
    falla "quedaron cambios sin guardar en Git"; }

echo
echo "================================================================"
echo " ANALISIS FINAL COMPLETO. Reporte: reporte/reporte.pdf"
echo " Historial:  git log --oneline | head"
echo "================================================================"
