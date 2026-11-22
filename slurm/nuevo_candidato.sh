#!/bin/bash
# =============================================================================
#  nuevo_candidato.sh - crea un candidato extra con otra malla de k-points,
#                       lo guarda en Git y lo envia. Solo se usa si el resumen
#                       de la calibracion lo recomienda.
#  Uso:  bash slurm/nuevo_candidato.sh K8 10 10 3
# =============================================================================
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
falla() { echo "  [FALLA] $1"; exit 1; }

[ $# -eq 4 ] || falla "uso: bash slurm/nuevo_candidato.sh NOMBRE N1 N2 N3   (ejemplo: K8 10 10 3)"
K="$1"; N1="$2"; N2="$3"; N3="$4"
echo "$K" | grep -qE '^K[0-9]+$' || falla "el nombre debe ser K seguido de un numero (K8, K9, ...)"
for n in "$N1" "$N2" "$N3"; do
    echo "$n" | grep -qE '^[1-9][0-9]*$' || falla "las divisiones de la malla deben ser enteros positivos"
done
grep -qE "^$K[[:space:]]" entradas/candidatos.txt && falla "$K ya existe en entradas/candidatos.txt"
[ -s entradas/POSCAR_bulk_2x2x1_VSe ] || falla "falta el POSCAR (corre primero: bash slurm/preparar.sh)"

PY="$(command -v python3 || command -v python)"
NK=$("$PY" analisis/resumen_t0.py --nkpts "$N1" "$N2" "$N3" 2>/dev/null)
[ -n "$NK" ] || NK="?"
printf "%-3s %10s %3s %3s %6s  %s\n" "$K" "$N1" "$N2" "$N3" "$NK" "$((2 * N1))x$((2 * N2))x$N3" >> entradas/candidatos.txt
D="entradas/candidatos/$K"
mkdir -p "$D"
cp entradas/INCAR_base "$D/INCAR"
cp entradas/POSCAR_bulk_2x2x1_VSe "$D/POSCAR"
printf "Malla Gamma %sx%sx%s de la supercelda 2x2x1 (= %sx%sx%s en la primitiva)\n0\nGamma\n%s %s %s\n0 0 0\n" \
       "$N1" "$N2" "$N3" "$((2 * N1))" "$((2 * N2))" "$N3" "$N1" "$N2" "$N3" > "$D/KPOINTS"
git add entradas/candidatos.txt "$D" && git commit -q -m "Tarea 0: candidato $K (malla ${N1}x${N2}x${N3})" \
    || falla "no se pudo hacer el commit"
echo "  [OK] $K creado: malla ${N1}x${N2}x${N3}, $NK k-points irreducibles esperados (commit $(git rev-parse --short HEAD))"
bash slurm/enviar.sh "$K"
