#!/bin/bash
# =============================================================================
#  parteA1_preparar.sh - verifica que todo este listo para la Tarea A.1
#                        (escalamiento fuerte con MPI) y guarda los scripts en Git.
#
#  Se ejecuta en el login, desde la carpeta del repositorio:
#      bash slurm/parteA1_preparar.sh
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
P_LISTA="1 2 4 8 12 16 24"

echo "Repositorio: $REPO"

# -----------------------------------------------------------------------------
paso "1/5  Archivos de la Parte A.1"
for f in slurm/vasp_job.sh slurm/parteA1_enviar.sh slurm/parteA1_estado.sh \
         analisis/parteA1_resumen.py README.md; do
    [ -s "$f" ] || falla "falta $f: el paquete no se descomprimio dentro de examen_silvera (ver la guia)"
done
grep -q "afinidad_ok" slurm/vasp_job.sh || falla "slurm/vasp_job.sh es la version anterior: descomprime de nuevo el paquete"
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
paso "2/5  Caso base de la Tarea 0"
for f in INCAR KPOINTS POSCAR ORIGEN.txt; do
    [ -s "entradas/caso_base/$f" ] || falla "falta entradas/caso_base/$f (la Tarea 0 no se fijo)"
done
leer() { awk -F= -v t="$1" '$1 ~ "^[ \t]*"t"[ \t]*$" {split($2,a,"#"); gsub(/[ \t]/,"",a[1]); print a[1]}' entradas/caso_base/INCAR; }
NCORE=$(leer NCORE); KPAR=$(leer KPAR); NBANDS=$(leer NBANDS); ENCUT=$(leer ENCUT)
MALLA=$(sed -n 4p entradas/caso_base/KPOINTS | tr -s ' ' | sed 's/^ //;s/ $//')
CAND=$(grep '^candidato' entradas/caso_base/ORIGEN.txt | cut -d= -f2 | tr -d ' ')
T1=$(grep '^t1_min' entradas/caso_base/ORIGEN.txt | cut -d= -f2 | tr -d ' ')
echo "       candidato $CAND, malla $(echo $MALLA | tr " " "x"), ENCUT = $ENCUT eV, NBANDS = $NBANDS, t1 (Tarea 0) = $T1 min"
[ "$NCORE" = "1" ] || falla "NCORE = $NCORE en el caso base (la A.1 usa NCORE = 1)"
[ "$KPAR" = "1" ]  || falla "KPAR = $KPAR en el caso base (la A.1 usa KPAR = 1)"
for p in $P_LISTA; do
    [ $((NBANDS % p)) -eq 0 ] || falla "NBANDS = $NBANDS no es divisible por p = $p"
done
ok "NCORE = 1 y KPAR = 1 (valores por defecto de VASP); NBANDS = $NBANDS divisible por todos los p"

# -----------------------------------------------------------------------------
paso "3/5  POTCAR y librerias"
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
if [ -s entradas/POTCAR.sha256 ]; then
    A=$(sha256sum entradas/POTCAR | awk '{print $1}')
    B=$(awk '{print $1}' entradas/POTCAR.sha256)
    [ "$A" = "$B" ] || falla "el POTCAR cambio desde la Tarea 0 (sha256 distinto)"
fi
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/ (corre bash slurm/preparar.sh)"
ok "POTCAR igual al de la Tarea 0 (mismo sha256) y librerias FFTW en lib/"

# -----------------------------------------------------------------------------
paso "4/5  Registro previo"
if [ -s logs/jobs_A1.txt ]; then
    aviso "ya existe logs/jobs_A1.txt: la A.1 ya se envio. Para ver el estado: bash slurm/parteA1_estado.sh"
else
    ok "la A.1 todavia no se ha enviado"
fi

# -----------------------------------------------------------------------------
paso "5/5  Git"
git add slurm analisis README.md || falla "git add fallo"
if git diff --cached --quiet; then
    ok "nada nuevo que guardar (los scripts ya estaban en Git)"
else
    git commit -q -m "Parte A.1: escalamiento fuerte con MPI (scripts y README)" || falla "git commit fallo"
    ok "commit $(git rev-parse --short HEAD)"
fi
[ -z "$(git status --porcelain -- entradas slurm analisis README.md)" ] || falla "quedaron cambios sin guardar en Git"

echo
echo "================================================================"
echo " PREPARACION A.1 COMPLETA. Siguiente paso:  bash slurm/parteA1_enviar.sh"
echo "================================================================"
