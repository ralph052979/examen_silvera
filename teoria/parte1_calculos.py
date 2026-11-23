#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parte1_calculos.py - calculos numericos de la Parte I (Problemas 1 y 3). El desarrollo
completo esta en teoria/parte1.md. Solo usa la biblioteca estandar; imprime ASCII.

Uso:  python3 teoria/parte1_calculos.py
"""
from __future__ import print_function, division


def amdahl(f, p):
    """aceleracion de Amdahl: S(p) = 1 / (f + (1 - f)/p)"""
    return 1.0 / (f + (1.0 - f) / p)


def gustafson(f, p):
    """aceleracion escalada de Gustafson: S(p) = f + (1 - f) p = p - f (p - 1)"""
    return f + (1.0 - f) * p


print("Problema 1 (b): f = 0.05")
f = 0.05
print("  S(24)  = 1/(0.05 + 0.95/24) = 1/%.7f = %.4f" % (f + (1 - f) / 24, amdahl(f, 24)))
print("  S(inf) = 1/f = %.1f" % (1 / f))
print("  E(24)  = S(24)/24 = %.4f" % (amdahl(f, 24) / 24))
print("  Comparacion (c): Gustafson con la misma fraccion serial: S(24) = 24 - 0.05*23 = %.2f" % gustafson(f, 24))
print()

print("Problema 3: t(1) = 100 s, t(4) = 30 s, t(16) = 11 s, t(24) = 9 s")
t = {1: 100.0, 4: 30.0, 16: 11.0, 24: 9.0}
print("  %-4s %8s %8s %8s %24s" % ("p", "t (s)", "S(p)", "E(p)", "Karp-Flatt e(p)"))
for p in sorted(t):
    S = t[1] / t[p]
    E = S / p
    kf = "-" if p == 1 else "%.4f" % ((1 / S - 1 / p) / (1 - 1 / p))
    print("  %-4d %8.1f %8.4f %8.4f %24s" % (p, t[p], S, E, kf))
# ajuste de minimos cuadrados de t(p) = a + b/p  (a = parte serial, b = parte paralelizable)
ps = sorted(t)
xs = [1.0 / p for p in ps]
ys = [t[p] for p in ps]
n = len(ps)
sx, sy = sum(xs), sum(ys)
sxx, sxy = sum(v * v for v in xs), sum(u * v for u, v in zip(xs, ys))
b = (n * sxy - sx * sy) / (n * sxx - sx * sx)
a = (sy - b * sx) / n
f_ajuste = a / (a + b)
print("  Ajuste t(p) = a + b/p (minimos cuadrados): a = %.3f s, b = %.3f s  ->  f = a/(a + b) = %.4f" % (a, b, f_ajuste))
# ajuste de Amdahl sobre 1/S - 1/p = f (1 - 1/p) (sin peso, p > 1)
X = [1 - 1.0 / p for p in ps if p > 1]
Y = [t[p] / t[1] - 1.0 / p for p in ps if p > 1]
f_lin = sum(u * v for u, v in zip(X, Y)) / sum(u * u for u in X)
print("  Ajuste lineal 1/S - 1/p = f (1 - 1/p) (p > 1):            f = %.4f" % f_lin)
for fe in (f_ajuste, f_lin):
    print("  Con f = %.4f, E(p) = 1/(f p + 1 - f) = 0.5 en p = (1 + f)/f = %.1f" % (fe, (1 + fe) / fe))
print("  Entre los p medidos, E < 0.5 por primera vez en p = 24 (E(16) = %.3f, E(24) = %.3f)" % (
    t[1] / t[16] / 16, t[1] / t[24] / 24))
