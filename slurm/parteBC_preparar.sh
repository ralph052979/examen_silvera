#!/bin/bash
# =============================================================================
#  parteBC_preparar.sh - Partes B y C.2: genera las entradas del escalamiento debil (B.1)
#                        y de la serie de ENCUT (C.2), verifica los scripts y lo guarda
#                        todo en Git. No usa VASP ni Slurm.
#
#  B.1 (escalamiento debil): el problema crece con p a traves de los k-points:
#  p procesos MPI con KPAR = p (cada proceso calcula sus propios 3 k-points, con las
#  144 bandas) y una malla Gamma con exactamente 3p k-points irreducibles. Las mallas
#  se eligieron contando los k-points irreducibles con la simetria de la supercelda
#  (el conteo reproduce los NKPTS que dio VASP para las 7 mallas de la Tarea 0):
#      p :  1      2      4      8      12      16       24
#    malla 3x3x1  3x3x2  5x5x3  6x6x4  10x10x3 14x14x2  10x10x6
#   NKPTS  3      6      12     24     36      48       72
#  El resto del INCAR es el del caso base (ENCUT = 400, NBANDS = 144, NCORE = 1).
#
#  C.2 (precision vs costo): ENCUT = 250 ... 800 eV con la malla 3x3x1 (3 k-points) en
#  1 nucleo; el resto del INCAR es el del caso base. Referencia: ENCUT = 800 eV.
#
#  B.2 (hibrido) usa entradas/caso_base/ sin cambios.
#
#  Uso (desde la carpeta del repositorio):   bash slurm/parteBC_preparar.sh
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
falla() { echo "  [FALLA] $1"; echo; echo "Mandame una captura."; exit 1; }
ok()    { echo "  [OK] $1"; }
paso()  { echo; echo "=== $1 ==="; }

MALLAS="1:3:3:1 2:3:3:2 4:5:5:3 8:6:6:4 12:10:10:3 16:14:14:2 24:10:10:6"
ENCUTS="250 300 350 400 450 500 550 600 700 800"

paso "1/5  Archivos de las Partes B y C"
for f in slurm/vasp_job.sh slurm/parteB_enviar.sh slurm/parteC_enviar.sh slurm/parteC2_serie.sh \
         slurm/parteBC_estado.sh analisis/parteB_resumen.py analisis/parteA3_graficas.py README.md; do
    [ -s "$f" ] || falla "falta $f: el paquete no se descomprimio dentro de examen_silvera (ver la guia)"
done
grep -q "mem_total_mb" slurm/vasp_job.sh && grep -q 'SUFIJO:-' slurm/vasp_job.sh \
    || falla "slurm/vasp_job.sh es una version anterior: descomprime de nuevo el paquete"
CR=$(grep -l $'\r' slurm/*.sh analisis/*.py README.md 2>/dev/null)
[ -z "$CR" ] || falla "estos archivos tienen finales de linea de Windows: $CR (no los abras en Windows)"
for f in slurm/*.sh; do bash -n "$f" || falla "error de sintaxis en $f"; done
for f in analisis/*.py; do python3 -m py_compile "$f" 2>/dev/null || falla "error de sintaxis en $f"; done
rm -rf analisis/__pycache__
ok "scripts completos y sin errores de sintaxis"

paso "2/5  Caso base, POTCAR y librerias"
for f in INCAR KPOINTS POSCAR; do
    [ -s "entradas/caso_base/$f" ] || falla "falta entradas/caso_base/$f"
done
[ -s entradas/POTCAR ] || falla "falta entradas/POTCAR"
[ -s lib/libfftw3_omp.so.3 ] && [ -s lib/libfftw3.so.3 ] || falla "faltan las librerias de lib/"
B=entradas/caso_base/INCAR
[ "$(grep -cE '^KPAR +=' "$B")" -eq 1 ] || falla "el INCAR del caso base no tiene exactamente una linea KPAR"
[ "$(grep -cE '^ENCUT +=' "$B")" -eq 1 ] || falla "el INCAR del caso base no tiene exactamente una linea ENCUT"
grep -qE '^NBANDS *= *144 *$' "$B" || falla "el INCAR del caso base no tiene NBANDS = 144"
grep -qE '^NCORE *= *1 *$' "$B" || falla "el INCAR del caso base no tiene NCORE = 1"
ok "caso base, POTCAR y librerias en orden"
if MPL=$(python3 -c "import matplotlib; print(matplotlib.__version__)" 2>/dev/null); then
    ok "matplotlib $MPL (para las graficas)"
else
    echo "  [AVISO] falta matplotlib (solo para las graficas; no impide enviar las corridas)."
    echo "          Instalalo una vez con:  python3 -m pip install --user matplotlib"
fi

paso "3/5  B.1: entradas del escalamiento debil (entradas/debil/P<p>/)"
for m in $MALLAS; do
    IFS=: read -r P N1 N2 N3 <<< "$m"
    D="entradas/debil/P$P"
    mkdir -p "$D"
    sed -e "s/^KPAR .*/KPAR   = $P/" \
        -e "s/^SYSTEM .*/SYSTEM = WSe2 bulk 2H 2x2x1 V_Se - Parte B.1 escalamiento debil p = $P/" "$B" > "$D/INCAR"
    printf "Malla Gamma %sx%sx%s de la supercelda 2x2x1 (escalamiento debil: %s k-points irreducibles para p = %s)\n0\nGamma\n%s %s %s\n0 0 0\n" \
           "$N1" "$N2" "$N3" "$((3 * P))" "$P" "$N1" "$N2" "$N3" > "$D/KPOINTS"
    cp entradas/caso_base/POSCAR "$D/POSCAR"
    [ "$(grep -c "^KPAR   = $P\$" "$D/INCAR")" -eq 1 ] || falla "no se pudo fijar KPAR = $P en $D/INCAR"
    [ "$(diff <(grep -vE '^(KPAR|SYSTEM) ' "$B") <(grep -vE '^(KPAR|SYSTEM) ' "$D/INCAR") | wc -l)" -eq 0 ] \
        || falla "$D/INCAR difiere del caso base en algo mas que KPAR"
    printf "  [OK] p = %-2s  malla %2s x %2s x %s  KPAR = %-2s -> %2s k-points irreducibles esperados (3 por proceso)\n" \
           "$P" "$N1" "$N2" "$N3" "$P" "$((3 * P))"
done
{
    echo "Escalamiento debil (Parte B.1): p procesos MPI, KPAR = p, 3 k-points irreducibles por proceso."
    echo "Generado por slurm/parteBC_preparar.sh a partir de entradas/caso_base/."
    echo "p malla nkpts_esperados"
    for m in $MALLAS; do IFS=: read -r P N1 N2 N3 <<< "$m"; echo "$P ${N1}x${N2}x${N3} $((3 * P))"; done
} > entradas/debil/ORIGEN.txt

paso "4/5  C.2: entradas de la serie de ENCUT (entradas/encut/E<ENCUT>/)"
for E in $ENCUTS; do
    D="entradas/encut/E$E"
    mkdir -p "$D"
    sed -e "s/^ENCUT .*/ENCUT  = $E/" \
        -e "s/^SYSTEM .*/SYSTEM = WSe2 bulk 2H 2x2x1 V_Se - Tarea C.2 ENCUT = $E eV, malla 3x3x1/" "$B" > "$D/INCAR"
    cp entradas/debil/P1/KPOINTS "$D/KPOINTS"
    cp entradas/caso_base/POSCAR "$D/POSCAR"
    [ "$(grep -c "^ENCUT  = $E\$" "$D/INCAR")" -eq 1 ] || falla "no se pudo fijar ENCUT = $E en $D/INCAR"
    [ "$(diff <(grep -vE '^(ENCUT|SYSTEM) ' "$B") <(grep -vE '^(ENCUT|SYSTEM) ' "$D/INCAR") | wc -l)" -eq 0 ] \
        || falla "$D/INCAR difiere del caso base en algo mas que ENCUT"
done
ok "ENCUT = $(echo $ENCUTS | sed 's/ /, /g') eV con la malla 3x3x1 (3 k-points); referencia 800 eV"
{
    echo "Serie de ENCUT (Tarea C.2): caso base con la malla Gamma 3x3x1 (3 k-points irreducibles), 1 nucleo."
    echo "Generado por slurm/parteBC_preparar.sh a partir de entradas/caso_base/. Referencia: ENCUT = 800 eV."
    echo "encut_eV: $ENCUTS"
} > entradas/encut/ORIGEN.txt

paso "5/5  Git"
git add entradas/debil entradas/encut slurm analisis README.md || falla "git add fallo"
if git diff --cached --quiet; then
    ok "nada nuevo que guardar (ya estaba en Git)"
else
    git commit -q -m "Partes B y C.2: entradas (KPAR = p con 3p k-points; serie de ENCUT), hibrido con memoria por proceso y graficas de la A.3" \
        || falla "git commit fallo"
    ok "commit $(git rev-parse --short HEAD)"
fi
[ -z "$(git status --porcelain -- entradas slurm analisis README.md)" ] || falla "quedaron cambios sin guardar en Git"

echo
echo "================================================================"
if [ -s logs/jobs_B1.txt ] || [ -s logs/jobs_B2.txt ]; then
    echo " PREPARACION COMPLETA. La Parte B ya se envio:  bash slurm/parteBC_estado.sh"
else
    echo " PREPARACION COMPLETA. Siguiente paso:  bash slurm/parteB_enviar.sh"
fi
echo "================================================================"
