#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generar_poscar_bulk.py - Tarea 0 del examen (Ruta A, VASP)

Construye la supercelda 2x2x1 de 2H-WSe2 en bulk y le quita UN atomo de Se
(vacancia V_Se). No relaja nada: es la geometria experimental.

  Estructura 2H (P6_3/mmc, No. 194), Schutte et al., J. Solid State Chem. 70, 207 (1987):
      a = 3.282 A, c = 12.96 A
      W  en 2c: (1/3, 2/3, 1/4), (2/3, 1/3, 3/4)
      Se en 4f: (1/3, 2/3, z), (2/3, 1/3, z+1/2), (2/3, 1/3, -z), (1/3, 2/3, 1/2-z), z = 0.621
      -> W-Se = 2.527 A y Se-Se (dentro de la capa) = 3.344 A

Escribe un POSCAR (formato VASP 5, orden de especies: W, Se) y un resumen
clave=valor con la geometria para la tabla de la Tarea 0. Verifica la estructura
y termina con codigo 1 si algo no cuadra.

Uso:
  python3 entradas/generar_poscar_bulk.py -o entradas/POSCAR_bulk_2x2x1_VSe \
                                          -r entradas/geometria_bulk_2x2x1_VSe.txt

Solo usa la biblioteca estandar de Python (no necesita numpy).
"""
from __future__ import print_function
import argparse
import math
import sys

A_EXP = 3.282     # A
C_EXP = 12.96     # A
Z_EXP = 0.621
N_SUPER = 2       # supercelda N x N x 1
ENLACE_MAX = 2.9  # A, corte para contar vecinos W-Se


def falla(msg):
    print("ERROR: " + msg)
    sys.exit(1)


def red_hexagonal(a, c):
    # convencion de 120 grados: a1 = a(1,0,0), a2 = a(-1/2, sqrt(3)/2, 0), a3 = (0,0,c)
    return [[a, 0.0, 0.0], [-a / 2.0, a * math.sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, c]]


def frac_a_cart(f, lat):
    return [f[0] * lat[0][k] + f[1] * lat[1][k] + f[2] * lat[2][k] for k in range(3)]


def dist_min_pbc(f1, f2, lat):
    """distancia minima entre dos posiciones fraccionarias con condiciones periodicas"""
    d0 = [f2[k] - f1[k] for k in range(3)]
    d0 = [x - math.floor(x + 0.5) for x in d0]
    mejor = 1e9
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                c = frac_a_cart([d0[0] + i, d0[1] + j, d0[2] + k], lat)
                mejor = min(mejor, math.sqrt(c[0] ** 2 + c[1] ** 2 + c[2] ** 2))
    return mejor


def envolver(x):
    """lleva x a [0, 1) sin dejar valores como 0.9999999999999999"""
    y = x - math.floor(x)
    if y > 1.0 - 1e-12:
        y = 0.0
    return y


def main():
    p = argparse.ArgumentParser(description="Supercelda 2x2x1 de 2H-WSe2 bulk con una vacancia de Se")
    p.add_argument("-o", "--salida", required=True, help="POSCAR de salida")
    p.add_argument("-r", "--resumen", default=None, help="archivo clave=valor con la geometria")
    p.add_argument("--a", type=float, default=A_EXP, help="parametro a (A), por defecto %.3f" % A_EXP)
    p.add_argument("--c", type=float, default=C_EXP, help="parametro c (A), por defecto %.2f" % C_EXP)
    p.add_argument("--z", type=float, default=Z_EXP, help="z del Se (4f), por defecto %.3f" % Z_EXP)
    arg = p.parse_args()

    a, c, z, n = arg.a, arg.c, arg.z, N_SUPER

    # ---------------------------------------------------------------- celda primitiva
    w_prim = [(1.0 / 3, 2.0 / 3, 0.25), (2.0 / 3, 1.0 / 3, 0.75)]
    se_prim = [(1.0 / 3, 2.0 / 3, z), (2.0 / 3, 1.0 / 3, z + 0.5),
               (2.0 / 3, 1.0 / 3, -z), (1.0 / 3, 2.0 / 3, 0.5 - z)]

    # ---------------------------------------------------------------- supercelda n x n x 1
    def repetir(lista):
        out = []
        for (x, y, zz) in lista:
            for i in range(n):
                for j in range(n):
                    out.append((envolver((x + i) / n), envolver((y + j) / n), envolver(zz)))
        return out

    w_sc = repetir(w_prim)
    se_sc = repetir(se_prim)
    lat = red_hexagonal(n * a, c)

    if len(w_sc) != 8 or len(se_sc) != 16:
        falla("la supercelda pristina deberia tener 8 W y 16 Se (tiene %d y %d)" % (len(w_sc), len(se_sc)))

    # ---------------------------------------------------------------- vacancia
    # Se quita el primer Se de la lista: (1/6, 1/3, z), el Se inferior de la segunda capa.
    # En el bulk pristino todos los Se son equivalentes por simetria.
    se_quitado = se_sc.pop(0)
    se_quitado_cart = frac_a_cart(se_quitado, lat)

    # ---------------------------------------------------------------- verificaciones
    atomos = [("W", f) for f in w_sc] + [("Se", f) for f in se_sc]
    nat = len(atomos)
    if nat != 23:
        falla("se esperaban 23 atomos y hay %d" % nat)

    dmin = 1e9
    for i in range(nat):
        for j in range(i + 1, nat):
            dmin = min(dmin, dist_min_pbc(atomos[i][1], atomos[j][1], lat))
    if dmin < 2.3:
        falla("hay atomos demasiado cerca (d_min = %.3f A)" % dmin)

    coord = []
    d_wse = []
    for fw in w_sc:
        vec = [dist_min_pbc(fw, fs, lat) for fs in se_sc]
        cerca = [d for d in vec if d < ENLACE_MAX]
        coord.append(len(cerca))
        d_wse.extend(cerca)
    n5 = sum(1 for k in coord if k == 5)
    n6 = sum(1 for k in coord if k == 6)
    if n5 != 3 or n6 != 5:
        falla("coordinacion inesperada: %d W con 5 Se y %d W con 6 Se (esperado 3 y 5)" % (n5, n6))

    vecinos_vac = sorted(dist_min_pbc(se_quitado, fw, lat) for fw in w_sc)[:3]
    # Se de la primera capa en z-1/2 y 1-z  ->  separacion vertical (3/2 - 2z) c
    se_se_intra = (1.5 - 2.0 * z) * c

    # ---------------------------------------------------------------- POSCAR
    lin = []
    lin.append("WSe2 2H bulk 2x2x1 V_Se  a=%.4f c=%.4f z=%.4f (exp. Schutte 1987) sin relajar" % (a, c, z))
    lin.append("   1.0")
    for v in lat:
        lin.append("  %22.16f %22.16f %22.16f" % (v[0], v[1], v[2]))
    lin.append("   W   Se")
    lin.append("   %d   %d" % (len(w_sc), len(se_sc)))
    lin.append("Direct")
    for nombre, f in atomos:
        lin.append("  %20.16f %20.16f %20.16f   %s" % (f[0], f[1], f[2], nombre))
    with open(arg.salida, "w") as fh:
        fh.write("\n".join(lin) + "\n")

    vol = abs(lat[0][0] * lat[1][1] * lat[2][2])
    res = [
        ("estructura", "2H-WSe2 bulk, P6_3/mmc, supercelda 2x2x1, 1 vacancia de Se"),
        ("fuente", "Schutte, de Boer y Jellinek, J. Solid State Chem. 70, 207 (1987)"),
        ("a_primitiva_A", "%.4f" % a),
        ("c_A", "%.4f" % c),
        ("z_Se", "%.4f" % z),
        ("a_supercelda_A", "%.4f" % (n * a)),
        ("volumen_A3", "%.3f" % vol),
        ("n_atomos", "%d" % nat),
        ("n_W", "%d" % len(w_sc)),
        ("n_Se", "%d" % len(se_sc)),
        ("d_min_A", "%.4f" % dmin),
        ("d_W_Se_A", "%.4f" % (sum(d_wse) / len(d_wse))),
        ("d_Se_Se_intracapa_A", "%.4f" % se_se_intra),
        ("W_con_5_Se", "%d" % n5),
        ("W_con_6_Se", "%d" % n6),
        ("Se_quitado_frac", "%.6f %.6f %.6f" % se_quitado),
        ("Se_quitado_cart_A", "%.4f %.4f %.4f" % tuple(se_quitado_cart)),
        ("W_vecinos_de_la_vacancia_A", " ".join("%.4f" % d for d in vecinos_vac)),
    ]
    if arg.resumen:
        with open(arg.resumen, "w") as fh:
            for k, v in res:
                fh.write("%s = %s\n" % (k, v))
    for k, v in res:
        print("  %-28s %s" % (k, v))
    print("OK: POSCAR escrito en " + arg.salida)


if __name__ == "__main__":
    main()
