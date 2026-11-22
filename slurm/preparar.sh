#!/bin/bash
# =============================================================================
#  preparar.sh - prepara el repositorio para la Tarea 0 (WSe2 bulk 2x2x1 con V_Se).
#
#  Se ejecuta en el login, desde la carpeta del repositorio:
#      bash slurm/preparar.sh
#  Solo copia y verifica archivos, genera el POSCAR (menos de 1 s) y hace un
#  commit. NO ejecuta VASP. Se puede repetir sin riesgo.
#
#  Variables opcionales:
#    POTCAR_ORIGEN  ruta a un POTCAR (W_sv y luego Se) si no esta en los sitios usuales
#    ANTERIOR       carpeta del intento anterior (por defecto ../examen_silvera_monocapa)
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
PADRE="$(dirname "$REPO")"
ANTERIOR="${ANTERIOR:-$PADRE/examen_silvera_monocapa}"
TESIS_POTCAR="${TESIS_POTCAR:-/work/physics/users/rsilvera/CALCULOS_WSe2_AC/04-Monocapa/2-Relax/2-estatico/POTCAR}"
LIBSRC="${LIBSRC:-/lib64}"

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
paso "1/8  Carpeta correcta"
for x in entradas/tesis_monocapa slurm/correr_vasp.sh slurm/tarea0_preparar.sh \
         analisis/tarea0_resumen.py analisis/parteA_analisis.py; do
    [ -e "$x" ] && falla "esta carpeta tiene archivos del intento anterior ($x). Ver el Paso 0 de la guia."
done
[ "$(basename "$REPO")" = "examen_silvera" ] || aviso "la carpeta se llama $(basename "$REPO"); el examen pide examen_<apellido>"
ok "carpeta limpia (solo archivos de esta version)"

# -----------------------------------------------------------------------------
paso "2/8  Archivos de texto con fin de linea Unix"
ARCH="slurm/*.sh analisis/*.py entradas/*.py entradas/INCAR_base entradas/candidatos.txt README.md .gitignore"
CR=$(grep -l $'\r' $ARCH 2>/dev/null)
if [ -n "$CR" ]; then
    sed -i 's/\r$//' $CR
    ok "corregido el fin de linea de Windows en: $(echo $CR)"
else
    ok "sin fin de linea de Windows (CRLF)"
fi
NOASCII=$(LC_ALL=C grep -l $'[\x80-\xFF]' entradas/INCAR_base entradas/candidatos.txt 2>/dev/null)
[ -z "$NOASCII" ] || falla "hay caracteres no ASCII (tildes) en: $NOASCII. VASP podria leerlos mal."

# -----------------------------------------------------------------------------
paso "3/8  Python"
PY="$(command -v python3 || command -v python || true)"
[ -n "$PY" ] || falla "no se encontro python3 ni python"
"$PY" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 5) else 1)' \
    || falla "se necesita Python >= 3.5 ($("$PY" --version 2>&1))"
ok "$("$PY" --version 2>&1)  ($PY)"

# -----------------------------------------------------------------------------
paso "4/8  POTCAR (debe ser W_sv y luego Se, PAW PBE)"
ORIG=""
for cand in "${POTCAR_ORIGEN:-}" "$REPO/entradas/POTCAR" "$ANTERIOR/entradas/POTCAR" "$TESIS_POTCAR"; do
    if [ -n "$cand" ] && [ -s "$cand" ]; then ORIG="$cand"; break; fi
done
[ -n "$ORIG" ] || falla "no encontre un POTCAR. Copia el de tu tesis (W_sv + Se) con: cp RUTA/POTCAR entradas/POTCAR  y repite este paso"
[ "$ORIG" = "$REPO/entradas/POTCAR" ] || cp "$ORIG" entradas/POTCAR || falla "no se pudo copiar $ORIG"
echo "       origen: $ORIG"
NT=$(grep -c "TITEL" entradas/POTCAR)
[ "$NT" -eq 2 ] || falla "el POTCAR tiene $NT bloques (se esperaban 2: W_sv y Se)"
T1=$(grep "TITEL" entradas/POTCAR | sed -n 1p | sed 's/.*= *//')
T2=$(grep "TITEL" entradas/POTCAR | sed -n 2p | sed 's/.*= *//')
echo "       bloque 1: $T1"
echo "       bloque 2: $T2"
echo "$T1" | grep -qE "(^| )W_sv( |$)" || falla "el primer bloque no es W_sv ($T1)"
echo "$T2" | grep -qE "(^| )Se( |$)"   || falla "el segundo bloque no es Se ($T2)"
ZV=$(awk '/TITEL/{n++} /ZVAL/{ if (!(n in z)) { s=$0; sub(/.*ZVAL *= */,"",s); split(s,a," "); z[n]=a[1] } } END{print z[1], z[2]}' entradas/POTCAR)
EM=$(awk '/TITEL/{n++} /ENMAX/{ if (!(n in e)) { s=$0; sub(/.*ENMAX *= */,"",s); split(s,a,";"); gsub(/ /,"",a[1]); e[n]=a[1] } } END{print e[1], e[2]}' entradas/POTCAR)
Z1=${ZV%% *}; Z2=${ZV##* }; E1=${EM%% *}; E2=${EM##* }
[ "${Z1%%.*}" = "14" ] || falla "ZVAL del W_sv = $Z1 (se esperaba 14)"
[ "${Z2%%.*}" = "6" ]  || falla "ZVAL del Se = $Z2 (se esperaba 6)"
NPE=$(grep -cE "LEXCH *= *PE" entradas/POTCAR)
[ "$NPE" -eq 2 ] || falla "LEXCH no es PE en los dos bloques: este POTCAR no es PBE (el INCAR usa GGA = PE)"
{
    echo "TITEL 1 = $T1"; echo "TITEL 2 = $T2"
    echo "ZVAL    = $Z1 $Z2"; echo "ENMAX   = $E1 $E2"; echo "LEXCH   = PE PE"
} > entradas/POTCAR.titel
sha256sum entradas/POTCAR | awk '{print $1"  POTCAR"}' > entradas/POTCAR.sha256
ok "PAW PBE, ZVAL = $Z1 y $Z2, ENMAX = $E1 y $E2 eV"
ok "identidad guardada en entradas/POTCAR.titel y POTCAR.sha256 (el POTCAR no va a Git)"

# -----------------------------------------------------------------------------
paso "5/8  POSCAR: supercelda 2x2x1 de WSe2 bulk con una vacancia de Se"
rm -f entradas/POSCAR_bulk_2x2x1_VSe entradas/geometria_bulk_2x2x1_VSe.txt
SALIDA=$("$PY" entradas/generar_poscar_bulk.py -o entradas/POSCAR_bulk_2x2x1_VSe \
         -r entradas/geometria_bulk_2x2x1_VSe.txt 2>&1)
RCG=$?
echo "$SALIDA" | sed 's/^/     /'
[ "$RCG" -eq 0 ] || falla "el generador del POSCAR fallo (mensaje arriba)"
[ -s entradas/POSCAR_bulk_2x2x1_VSe ] || falla "no se genero el POSCAR (mensaje arriba)"
L6=$(sed -n 6p entradas/POSCAR_bulk_2x2x1_VSe | tr -s ' ' | sed 's/^ //;s/ $//')
L7=$(sed -n 7p entradas/POSCAR_bulk_2x2x1_VSe | tr -s ' ' | sed 's/^ //;s/ $//')
[ "$L6" = "W Se" ] || falla "linea 6 del POSCAR = '$L6' (se esperaba 'W Se', mismo orden que el POTCAR)"
[ "$L7" = "8 15" ] || falla "linea 7 del POSCAR = '$L7' (se esperaba '8 15')"
grep -q "^W_con_5_Se = 3$" entradas/geometria_bulk_2x2x1_VSe.txt || falla "la vacancia no quedo bien (deben ser 3 W con 5 Se)"
ok "23 atomos (8 W + 15 Se), 3 W vecinos de la vacancia, orden W Se igual al POTCAR"

# -----------------------------------------------------------------------------
paso "6/8  INCAR y candidatos de calibracion"
ENCUT=$(awk -F= '/^[ \t]*ENCUT[ \t]*=/{split($2,a,"#"); gsub(/[ \t]/,"",a[1]); print a[1]}' entradas/INCAR_base)
NBANDS=$(awk -F= '/^[ \t]*NBANDS[ \t]*=/{split($2,a,"#"); gsub(/[ \t]/,"",a[1]); print a[1]}' entradas/INCAR_base)
[ -n "$ENCUT" ] && [ -n "$NBANDS" ] || falla "no pude leer ENCUT o NBANDS de entradas/INCAR_base"
EMAX=$(awk -v a="$E1" -v b="$E2" 'BEGIN{print (a>b)?a:b}')
awk -v e="$ENCUT" -v m="$EMAX" 'BEGIN{exit !(e+0 >= m+0)}' || falla "ENCUT = $ENCUT es menor que el ENMAX del POTCAR ($EMAX)"
NELECT=$(awk -v a="$Z1" -v b="$Z2" 'BEGIN{printf "%d", 8*a + 15*b}')
[ $((NBANDS % 48)) -eq 0 ] || falla "NBANDS = $NBANDS no es multiplo de 48"
[ "$NBANDS" -ge $((NELECT / 2 + 20)) ] || falla "NBANDS = $NBANDS es muy chico para NELECT = $NELECT"
ok "ENCUT = $ENCUT eV (ENMAX maximo = $EMAX), NELECT = $NELECT, NBANDS = $NBANDS (multiplo de 48)"
rm -rf entradas/candidatos
NCAND=0
while read -r K N1 N2 N3 NK EQ; do
    case "$K" in ''|\#*) continue ;; esac
    D="entradas/candidatos/$K"
    mkdir -p "$D"
    cp entradas/INCAR_base "$D/INCAR"
    cp entradas/POSCAR_bulk_2x2x1_VSe "$D/POSCAR"
    printf "Malla Gamma %sx%sx%s de la supercelda 2x2x1 (= %s en la primitiva)\n0\nGamma\n%s %s %s\n0 0 0\n" \
           "$N1" "$N2" "$N3" "$EQ" "$N1" "$N2" "$N3" > "$D/KPOINTS"
    printf "  [OK] %-3s malla %2s x %2s x %s  -> %3s k-points irreducibles esperados  (= %s en la primitiva)\n" \
           "$K" "$N1" "$N2" "$N3" "$NK" "$EQ"
    NCAND=$((NCAND + 1))
done < entradas/candidatos.txt
[ "$NCAND" -ge 1 ] || falla "entradas/candidatos.txt no tiene candidatos"

# -----------------------------------------------------------------------------
paso "7/8  Librerias FFTW que faltan en los nodos de computo"
mkdir -p lib logs corridas
for L in libfftw3_omp.so.3 libfftw3.so.3; do
    if [ -e "$LIBSRC/$L" ]; then
        cp -L "$LIBSRC/$L" lib/ || falla "no se pudo copiar $LIBSRC/$L"
    elif [ -e "$ANTERIOR/lib/$L" ]; then
        cp -L "$ANTERIOR/lib/$L" lib/ || falla "no se pudo copiar $ANTERIOR/lib/$L"
    else
        falla "no encontre $L en $LIBSRC ni en $ANTERIOR/lib"
    fi
done
ok "copiadas a lib/ (cada job las agrega a LD_LIBRARY_PATH)"

# -----------------------------------------------------------------------------
paso "8/8  Git"
command -v git >/dev/null || falla "git no esta disponible"
if [ ! -d .git ]; then
    git init -q . || falla "git init fallo"
    git symbolic-ref HEAD refs/heads/main
    ok "repositorio creado (rama main)"
fi
git config user.name >/dev/null || git config user.name "Rafael Silvera Altamiranda"
if ! git config user.email >/dev/null; then
    EMAIL=$(git -C "$ANTERIOR" config user.email 2>/dev/null)
    if [ -z "$EMAIL" ]; then
        read -r -p "  Escribe tu correo institucional para los commits: " EMAIL
    fi
    [ -n "$EMAIL" ] || falla "se necesita un correo para git"
    git config user.email "$EMAIL"
fi
ok "autor: $(git config user.name) <$(git config user.email)>"
git add -A . || falla "git add fallo"
git ls-files --cached | grep -qE "(^|/)POTCAR$" && falla "el POTCAR quedo en Git; revisa .gitignore"
git ls-files --cached | grep -q "^lib/" && falla "lib/ quedo en Git; revisa .gitignore"
if git diff --cached --quiet; then
    ok "nada nuevo que guardar (ya estaba todo en Git)"
else
    git commit -q -m "Tarea 0: WSe2 bulk 2x2x1 con V_Se (SCF), candidatos de calibracion y scripts" \
        || falla "git commit fallo"
    ok "commit $(git rev-parse --short HEAD)"
fi
echo "       fuera de Git (correcto): POTCAR, lib/, logs/, corridas/"

echo
echo "================================================================"
echo " PREPARACION COMPLETA. Siguiente paso:  bash slurm/enviar.sh"
echo "================================================================"
