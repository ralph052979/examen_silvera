#!/bin/bash
# =============================================================================
#  parteA2_preparar.sh - verifica que todo este listo para la Tarea A.2
#                        (OpenMP) y guarda los scripts nuevos en Git.
#
#  Se ejecuta en el login, desde la carpeta del repositorio:
#      bash slurm/parteA2_preparar.sh
#  No ejecuta VASP. Se puede repetir sin riesgo.
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1

paso()  { echo; echo "=== $1 ==="; }
ok()    { echo "  [OK] $1"; }
aviso() { echo "  [AVISO] $1"; }
falla() {
    echo "  [FALLA] $1"
    echo
    echo "La preparacion se detuvo aqui. No se envio nada. Mandame una captura de esta pantalla."
    exit 1
}

echo "Repositorio: $REPO"

# -----------------------------------------------------------------------------
paso "1/5  Archivos de la Parte A.2"
for f in slurm/vasp_job.sh slurm/parteA2_diagnostico.sh slurm/parteA2_enviar.sh slurm/parteA2_estado.sh \
         analisis/parteA2_diagnostico.py analisis/parteA2_resumen.py README.md; do
    [ -s "$f" ] || falla "falta $f: el paquete no se descomprimio dentro de examen_silvera (ver la guia)"
done
grep -q "hilos_muestras" slurm/vasp_job.sh || falla "slurm/vasp_job.sh es una version anterior: descomprime de nuevo el paquete"
CR=$(grep -l $'\r' slurm/*.sh analisis/*.py README.md 2>/dev/null)
if [ -n "$CR" ]; then
    sed -i 's/\r$//' $CR
    ok "corregido el fin de linea de Windows en: $(echo $CR)"
fi
for f in slurm/*.sh; do bash -n "$f" || falla "error de sintaxis en $f"; done
PY="$(command -v python3 || command -v python || true)"
[ -n "$PY" ] || falla "no se encontro python3"
for f in analisis/*.py; do
    "$PY" -c "import ast, sys; ast.parse(open(sys.argv[1]).read())" "$f" || falla "error de sintaxis en $f"
done
ok "scripts completos y sin errores de sintaxis"

# -----------------------------------------------------------------------------
paso "2/5  Caso base"
for f in INCAR KPOINTS POSCAR ORIGEN.txt; do
    [ -s "entradas/caso_base/$f" ] || falla "falta entradas/caso_base/$f"
done
leer() { awk -F= -v t="$1" '$1 ~ "^[ \t]*"t"[ \t]*$" {split($2,a,"#"); gsub(/[ \t]/,"",a[1]); print a[1]}' entradas/caso_base/INCAR; }
[ "$(leer NCORE)" = "1" ] || falla "NCORE no es 1 en el caso base (VASP lo exige con OpenMP)"
[ "$(leer KPAR)" = "1" ]  || falla "KPAR no es 1 en el caso base"
ok "caso base K7 sin cambios (NCORE = 1, KPAR = 1)"

# -----------------------------------------------------------------------------
paso "3/5  t(1) de la A.1"
[ -s logs/jobs_A1.txt ] || falla "no existe logs/jobs_A1.txt: la A.2 usa las corridas con p = 1 de la A.1"
N1=0
while read -r P R ID F H; do
    [ "$P" = "1" ] || continue
    E=$(grep '^estado=' "corridas/ex_A1_p1_r${R}_${ID}/tiempo.txt" 2>/dev/null | tail -1 | cut -d= -f2)
    [ "$E" = "OK" ] && N1=$((N1 + 1))
done < logs/jobs_A1.txt
[ "$N1" -ge 1 ] || falla "no hay corridas con p = 1 terminadas en la A.1"
ok "$N1 corridas con p = 1 de la A.1 disponibles como t(1)"

# -----------------------------------------------------------------------------
paso "4/5  POTCAR y librerias"
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
if [ -s entradas/POTCAR.sha256 ]; then
    A=$(sha256sum entradas/POTCAR | awk '{print $1}')
    B=$(awk '{print $1}' entradas/POTCAR.sha256)
    [ "$A" = "$B" ] || falla "el POTCAR cambio desde la Tarea 0 (sha256 distinto)"
fi
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
ok "POTCAR igual al de la Tarea 0 y librerias FFTW en lib/ (libfftw3_omp es la que usa OpenMP)"

# -----------------------------------------------------------------------------
paso "5/5  Git"
git add slurm analisis README.md || falla "git add fallo"
if git diff --cached --quiet; then
    ok "nada nuevo que guardar (los scripts ya estaban en Git)"
else
    git commit -q -m "Parte A.2: escalamiento fuerte con OpenMP (diagnostico, barrido, analisis)" || falla "git commit fallo"
    ok "commit $(git rev-parse --short HEAD)"
fi
[ -z "$(git status --porcelain -- entradas slurm analisis README.md)" ] || falla "quedaron cambios sin guardar en Git"

echo
echo "================================================================"
if [ -s logs/jobs_A2d.txt ]; then
    echo " PREPARACION A.2 COMPLETA. El diagnostico ya se envio:"
    echo "   python3 analisis/parteA2_diagnostico.py"
else
    echo " PREPARACION A.2 COMPLETA. Siguiente paso:  bash slurm/parteA2_diagnostico.sh"
fi
echo "================================================================"
