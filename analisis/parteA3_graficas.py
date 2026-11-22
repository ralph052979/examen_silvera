#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteA3_graficas.py - Tarea A.3: aceleracion S(p) y eficiencia E(p) de MPI (A.1) y de
OpenMP (A.2), con la curva ideal, la ley de Amdahl ajustada y un modelo con sobrecosto.

Modelos (p = nucleos):
  Amdahl          1/S = f + (1 - f)/p                    ->  1/S - 1/p = f (1 - 1/p)
  con sobrecosto  1/S = f + (1 - f)/p + c (p - 1)        ->  1/S - 1/p = f (1 - 1/p) + c (p - 1)
Ambos son lineales en sus parametros: se ajustan por minimos cuadrados ponderados con
la incertidumbre de S (sigma_y = sigma_S / S^2); si chi^2 reducido > 1 se agranda la
incertidumbre de los parametros por sqrt(chi^2 reducido). Con c > 0, S es maxima en
p* = sqrt((1 - f)/c). Para cada punto se da tambien la fraccion serial de Karp-Flatt
e(p) = (1/S - 1/p)/(1 - 1/p): si fuera solo la fraccion serial de Amdahl, e(p) seria
constante; si crece con p, el limite es un sobrecosto que crece con p (comunicacion,
sincronizacion, ancho de banda de memoria).

Se ajustan cuatro series: MPI y OpenMP con el tiempo de la corrida (A.1: SCF completo;
A.2: 4 iteraciones SCF) y MPI y OpenMP por iteracion SCF (LOOP del OUTCAR, primeras 4
iteraciones de ambas).

Lee:     datos/parteA1_resumen.csv, datos/parteA2_resumen.csv, datos/parteA2_mpi_vs_omp.csv
Escribe: reporte/figuras/A3_escalamiento.pdf y .png    S(p) y E(p), tiempo de la corrida
         reporte/figuras/A3_por_iteracion.pdf y .png   S(p) por iteracion y e(p) de Karp-Flatt
         datos/parteA3_ajustes.csv y .md               parametros y calidad de cada ajuste
         datos/parteA3_karp_flatt.csv                  e(p) de cada punto

Uso (desde la carpeta del repositorio):
  python3 analisis/parteA3_graficas.py             # ajustes, tablas y figuras
  python3 analisis/parteA3_graficas.py --guardar   # ademas hace commit de figuras y tablas

Necesita matplotlib solo para las figuras (los ajustes y las tablas se escriben igual).
Todo lo que imprime es ASCII.
"""
from __future__ import print_function, division
import argparse
import csv
import io
import math
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(REPO, "datos")
FIGS = os.path.join(REPO, "reporte", "figuras")
P_LISTA = [1, 2, 4, 8, 12, 16, 24]

# Colores (paleta categorica validada: azul y naranja; tintas neutras para texto y ejes)
C_MPI, C_OMP = "#2a78d6", "#eb6834"
TINTA, TINTA2, MUTED, REJILLA, EJE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"


# ----------------------------------------------------------------------------- lectura
def leer_csv(nombre):
    ruta = os.path.join(DATOS, nombre)
    if not os.path.isfile(ruta):
        return None
    with io.open(ruta, "r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def x(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def serie(filas, cs, cd, ct=None):
    """lista de (p, S, sigma_S, t) con los puntos que tienen S"""
    out = []
    for f in filas or []:
        p, s, sd = x(f.get("p")), x(f.get(cs)), x(f.get(cd))
        if p is None or s is None or s <= 0:
            continue
        out.append((int(p), s, sd if sd and sd > 0 else None, x(f.get(ct)) if ct else None))
    return sorted(out)


# ----------------------------------------------------------------------------- ajustes
def resolver2(a11, a12, a22, b1, b2):
    det = a11 * a22 - a12 * a12
    if abs(det) < 1e-300:
        return None
    return ((a22 * b1 - a12 * b2) / det, (a11 * b2 - a12 * b1) / det, a22 / det, -a12 / det, a11 / det)


def ajustar(pts):
    """pts: (p, S, sigma_S, t). Devuelve los dos ajustes y e(p) de Karp-Flatt"""
    q = [(p, s, sd) for p, s, sd, t in pts if p > 1]
    if len(q) < 2:
        return None
    sy_pos = sorted(sd / s ** 2 for p, s, sd in q if sd)
    piso = sy_pos[len(sy_pos) // 2] if sy_pos else 0.01
    X, Y, Z, W = [], [], [], []
    for p, s, sd in q:
        X.append(1.0 - 1.0 / p)
        Y.append(1.0 / s - 1.0 / p)
        Z.append(p - 1.0)
        sy = max(sd / s ** 2 if sd else piso, 0.25 * piso, 1e-6)
        W.append(1.0 / sy ** 2)
    n = len(X)
    r = {"n_puntos": n}
    # Amdahl: y = f x
    sxx = sum(w * a * a for w, a in zip(W, X))
    sxy = sum(w * a * b for w, a, b in zip(W, X, Y))
    f = sxy / sxx
    chi2 = sum(w * (b - f * a) ** 2 for w, a, b in zip(W, X, Y))
    r["chi2r_amdahl"] = chi2 / (n - 1) if n > 1 else float("nan")
    r["f"] = f
    r["f_desv"] = math.sqrt(1.0 / sxx) * max(1.0, math.sqrt(r["chi2r_amdahl"]))
    # Amdahl + sobrecosto: y = f x + c z
    sol = resolver2(sxx, sum(w * a * c for w, a, c in zip(W, X, Z)), sum(w * c * c for w, c in zip(W, Z)),
                    sxy, sum(w * c * b for w, c, b in zip(W, Z, Y)))
    if sol and n > 2:
        f2, c2, v11, v12, v22 = sol
        chi2b = sum(w * (b - f2 * a - c2 * cz) ** 2 for w, a, b, cz in zip(W, X, Y, Z))
        cr = chi2b / (n - 2)
        k = max(1.0, cr)
        r.update({"f2": f2, "f2_desv": math.sqrt(v11 * k), "c": c2, "c_desv": math.sqrt(v22 * k), "chi2r_sobrecosto": cr})
        if c2 > 0 and f2 < 1:
            ps = math.sqrt((1.0 - f2) / c2)
            r["p_max"] = ps
            r["S_max"] = 1.0 / (f2 + (1.0 - f2) / ps + c2 * (ps - 1.0))
    r["karp_flatt"] = [(p, (1.0 / s - 1.0 / p) / (1.0 - 1.0 / p),
                        (sd / s ** 2) / (1.0 - 1.0 / p) if sd else None) for p, s, sd in q]
    return r


def s_amdahl(p, f):
    return 1.0 / (f + (1.0 - f) / p)


def s_sobrecosto(p, f, c):
    return 1.0 / (f + (1.0 - f) / p + c * (p - 1.0))


# ----------------------------------------------------------------------------- figuras
def estilo(plt):
    plt.rcParams.update({
        "font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "legend.fontsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.edgecolor": EJE, "axes.linewidth": 0.8,
        "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2, "text.color": TINTA,
        "grid.color": REJILLA, "grid.linewidth": 0.6, "axes.grid": True, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 200, "figure.dpi": 100})


def eje_p(ax):
    try:
        ax.set_xscale("log", base=2)
    except (TypeError, ValueError):          # matplotlib < 3.3
        ax.set_xscale("log", basex=2)
    ax.set_xticks(P_LISTA)
    ax.set_xticklabels([str(p) for p in P_LISTA])
    ax.minorticks_off()
    ax.set_xlim(0.85, 30)
    ax.set_xlabel(u"p (núcleos)")


def eje_log_y(ax, ticks, lim):
    ax.set_yscale("log")
    ax.set_yticks(ticks)
    ax.set_yticklabels([("%g" % t) for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(*lim)


def puntos(ax, pts, color, marca, etiqueta, relleno=True, eficiencia=False):
    ps = [p for p, s, sd, t in pts]
    ys = [(s / p if eficiencia else s) for p, s, sd, t in pts]
    es = [((sd or 0.0) / p if eficiencia else (sd or 0.0)) for p, s, sd, t in pts]
    ax.errorbar(ps, ys, yerr=es, fmt=marca, ms=5.5, color=color, mfc=color if relleno else "white",
                mec=color, mew=1.3, elinewidth=0.9, capsize=2.5, label=etiqueta, zorder=4)


def curvas(ax, fit, color, modo, eficiencia=False, sobrecosto=True):
    if not fit:
        return
    pp = [1.0 + i * 0.05 for i in range(int(23.0 / 0.05) + 1)]
    ya = [s_amdahl(p, fit["f"]) / (p if eficiencia else 1.0) for p in pp]
    ax.plot(pp, ya, "-", color=color, lw=1.4, alpha=0.9, zorder=3,
            label=u"%s, Amdahl: f = %.3f ± %.3f" % (modo, fit["f"], fit["f_desv"]))
    if sobrecosto and "c" in fit:
        yb = [s_sobrecosto(p, fit["f2"], fit["c"]) / (p if eficiencia else 1.0) for p in pp]
        ax.plot(pp, yb, ":", color=color, lw=1.6, zorder=3,
                label=u"%s, con sobrecosto: f = %.3f, c = %.1e" % (modo, fit["f2"], fit["c"]))


def ideal(ax, eficiencia=False):
    pp = [1.0, 24.0]
    ax.plot(pp, [1.0, 1.0] if eficiencia else pp, "--", color=MUTED, lw=1.0, zorder=2, label="ideal")


def rotulo(ax, pts, texto, color_txt, eficiencia=False, dy=1.0):
    if not pts:
        return
    p, s, sd, t = pts[-1]
    y = s / p if eficiencia else s
    ax.annotate(texto, (p, y), xytext=(6, 0), textcoords="offset points", va="center", ha="left",
                fontsize=8, color=color_txt)


def letra(ax, t):
    ax.set_title(t, loc="left", fontsize=9, fontweight="bold", color=TINTA)


def ordenar(h, l, orden):
    """reordena las entradas de la leyenda segun el comienzo de cada etiqueta"""
    par = []
    for pref in orden:
        par += [(a, b) for a, b in zip(h, l) if b.startswith(pref) and (a, b) not in par]
    par += [(a, b) for a, b in zip(h, l) if (a, b) not in par]
    return [a for a, b in par], [b for a, b in par]


def figuras(plt, mpi, omp, fmpi, fomp, mpi_i, omp_i, fmpi_i, fomp_i):
    hechas = []
    os.path.isdir(FIGS) or os.makedirs(FIGS)
    # ---- figura 1: S(p) y E(p) con el tiempo de la corrida
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.5))
    ideal(a1)
    curvas(a1, fmpi, C_MPI, "MPI")
    curvas(a1, fomp, C_OMP, "OpenMP")
    puntos(a1, mpi, C_MPI, "o", u"MPI (A.1), medido")
    puntos(a1, omp, C_OMP, "s", u"OpenMP (A.2), medido")
    rotulo(a1, mpi, "MPI", TINTA2)
    rotulo(a1, omp, "OpenMP", TINTA2)
    eje_p(a1)
    eje_log_y(a1, [0.5, 1, 2, 4, 8, 16, 24], (0.5, 30))
    a1.set_ylabel(u"aceleración S(p) = t(1)/t(p)")
    letra(a1, "(a)")
    ideal(a2, True)
    curvas(a2, fmpi, C_MPI, "MPI", eficiencia=True)
    curvas(a2, fomp, C_OMP, "OpenMP", eficiencia=True)
    puntos(a2, mpi, C_MPI, "o", u"MPI (A.1), medido", eficiencia=True)
    puntos(a2, omp, C_OMP, "s", u"OpenMP (A.2), medido", eficiencia=True)
    rotulo(a2, mpi, "MPI", TINTA2, True)
    rotulo(a2, omp, "OpenMP", TINTA2, True)
    eje_p(a2)
    a2.set_ylim(0, 1.12)
    a2.set_ylabel(u"eficiencia E(p) = S(p)/p")
    letra(a2, "(b)")
    h, l = ordenar(*(a1.get_legend_handles_labels() + (["MPI (A.1)", "MPI, A", "MPI, c", "ideal",
                                                        "OpenMP (A.2)", "OpenMP, A", "OpenMP, c"],)))
    fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.01), handlelength=2.6, columnspacing=1.5)
    fig.tight_layout(rect=(0, 0.21, 1, 1))
    for ext in ("pdf", "png"):
        ruta = os.path.join(FIGS, "A3_escalamiento." + ext)
        fig.savefig(ruta, bbox_inches="tight")
        hechas.append(ruta)
    plt.close(fig)
    # ---- figura 2: S(p) por iteracion SCF y fraccion serial de Karp-Flatt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.6))
    ideal(a1)
    curvas(a1, fmpi_i, C_MPI, "MPI", sobrecosto=False)
    curvas(a1, fomp_i, C_OMP, "OpenMP", sobrecosto=False)
    puntos(a1, mpi_i, C_MPI, "o", u"MPI, por iteración SCF", relleno=False)
    puntos(a1, omp_i, C_OMP, "s", u"OpenMP, por iteración SCF", relleno=False)
    rotulo(a1, mpi_i, "MPI", TINTA2)
    rotulo(a1, omp_i, "OpenMP", TINTA2)
    eje_p(a1)
    eje_log_y(a1, [0.5, 1, 2, 4, 8, 16, 24], (0.5, 30))
    a1.set_ylabel(u"aceleración por iteración SCF")
    letra(a1, "(a)")
    for fit, color, marca, relleno, et in ((fmpi, C_MPI, "o", True, u"MPI, corrida"),
                                          (fmpi_i, C_MPI, "o", False, u"MPI, por iteración"),
                                          (fomp, C_OMP, "s", True, u"OpenMP, corrida"),
                                          (fomp_i, C_OMP, "s", False, u"OpenMP, por iteración")):
        if not fit:
            continue
        kf = [(p, e, se) for p, e, se in fit["karp_flatt"] if e > 0]
        if not kf:
            continue
        a2.errorbar([p for p, e, se in kf], [e for p, e, se in kf], yerr=[se or 0.0 for p, e, se in kf],
                    fmt=marca + "-", ms=5, lw=1.0, color=color, mfc=color if relleno else "white", mec=color, mew=1.2,
                    elinewidth=0.8, capsize=2, label=et, zorder=4)
    eje_p(a2)
    a2.set_yscale("log")
    todos = [e for fit in (fmpi, fmpi_i, fomp, fomp_i) if fit for p, e, se in fit["karp_flatt"] if e > 0]
    if todos:
        marcas = [t for t in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5) if min(todos) / 1.6 <= t <= max(todos) * 1.6]
        a2.set_yticks(marcas)
        a2.set_yticklabels(["%g" % t for t in marcas])
        a2.set_ylim(min(todos) / 1.8, max(todos) * 1.8)
    a2.minorticks_off()
    a2.set_ylabel(u"fracción serial de Karp-Flatt e(p)")
    letra(a2, "(b)")
    h1, l1 = ordenar(*(a1.get_legend_handles_labels() + (["MPI, por", "MPI, A", "ideal", "OpenMP, por", "OpenMP, A"],)))
    a1.legend(h1, l1, loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=2, handlelength=2.2, columnspacing=1.0,
              fontsize=7)
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=2, handlelength=2.2, columnspacing=1.0, fontsize=7)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        ruta = os.path.join(FIGS, "A3_por_iteracion." + ext)
        fig.savefig(ruta, bbox_inches="tight")
        hechas.append(ruta)
    plt.close(fig)
    return hechas


# ----------------------------------------------------------------------------- tablas
def tablas(series):
    cols = ["serie", "n_puntos", "f", "f_desv", "chi2r_amdahl", "S24_amdahl", "S_inf_amdahl",
            "f2", "f2_desv", "c", "c_desv", "chi2r_sobrecosto", "p_max", "S_max"]
    with io.open(os.path.join(DATOS, "parteA3_ajustes.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for nombre, fit in series:
            if not fit:
                continue
            d = dict(fit)
            d["serie"] = nombre
            d["S24_amdahl"] = s_amdahl(24.0, fit["f"])
            d["S_inf_amdahl"] = (1.0 / fit["f"]) if fit["f"] > 0 else None
            w.writerow([("" if d.get(k) is None else (("%.6g" % d[k]) if isinstance(d[k], float) else d[k])) for k in cols])
    with io.open(os.path.join(DATOS, "parteA3_karp_flatt.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["serie", "p", "e", "e_desv"])
        for nombre, fit in series:
            for p, e, se in (fit or {}).get("karp_flatt", []):
                w.writerow([nombre, p, "%.6g" % e, "" if se is None else "%.3g" % se])
    md = [u"# Tarea A.3 - ajustes de Amdahl", u"",
          u"Amdahl: 1/S = f + (1-f)/p. Con sobrecosto: 1/S = f + (1-f)/p + c(p-1). Minimos cuadrados ponderados "
          u"sobre 1/S - 1/p (p > 1); incertidumbres de 1 sigma, agrandadas por sqrt(chi2 reducido) si este pasa de 1.", u"",
          u"| serie | f (Amdahl) | chi2 red. | S(24) Amdahl | S max = 1/f | f (con sobrecosto) | c | chi2 red. | p* | S(p*) |",
          u"|---|---|---|---|---|---|---|---|---|---|"]
    for nombre, fit in series:
        if not fit:
            continue
        md.append(u"| %s | %.4f +/- %.4f | %.2f | %.2f | %s | %s | %s | %s | %s | %s |" % (
            nombre, fit["f"], fit["f_desv"], fit["chi2r_amdahl"], s_amdahl(24.0, fit["f"]),
            ("%.1f" % (1.0 / fit["f"])) if fit["f"] > 0 else "-",
            ("%.4f +/- %.4f" % (fit["f2"], fit["f2_desv"])) if "f2" in fit else "-",
            ("%.2e +/- %.1e" % (fit["c"], fit["c_desv"])) if "c" in fit else "-",
            ("%.2f" % fit["chi2r_sobrecosto"]) if "chi2r_sobrecosto" in fit else "-",
            ("%.1f" % fit["p_max"]) if "p_max" in fit else "-",
            ("%.2f" % fit["S_max"]) if "S_max" in fit else "-"))
    md += [u"", u"Fraccion serial de Karp-Flatt e(p) = (1/S - 1/p)/(1 - 1/p):", u"",
           u"| serie | " + u" | ".join(u"p = %d" % p for p in P_LISTA[1:]) + u" |",
           u"|---|" + u"---|" * (len(P_LISTA) - 1)]
    for nombre, fit in series:
        if not fit:
            continue
        kf = dict((p, e) for p, e, se in fit["karp_flatt"])
        md.append(u"| %s | " % nombre + u" | ".join(("%.3f" % kf[p]) if p in kf else u"-" for p in P_LISTA[1:]) + u" |")
    md.append(u"")
    with io.open(os.path.join(DATOS, "parteA3_ajustes.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md))


def imprimir(series):
    print("Tarea A.3: ajustes de Amdahl (p > 1, minimos cuadrados ponderados)")
    print()
    print("%-22s %3s %18s %7s %8s %18s %20s %7s %6s %6s" % ("serie", "n", "f Amdahl", "chi2r", "S(24)", "f sobrecosto",
                                                           "c", "chi2r", "p*", "S(p*)"))
    for nombre, fit in series:
        if not fit:
            print("%-22s   sin datos suficientes" % nombre)
            continue
        print("%-22s %3d %18s %7.2f %8.2f %18s %20s %7s %6s %6s" % (
            nombre, fit["n_puntos"], "%.4f +/- %.4f" % (fit["f"], fit["f_desv"]), fit["chi2r_amdahl"],
            s_amdahl(24.0, fit["f"]),
            ("%.4f +/- %.4f" % (fit["f2"], fit["f2_desv"])) if "f2" in fit else "-",
            ("%.2e +/- %.1e" % (fit["c"], fit["c_desv"])) if "c" in fit else "-",
            ("%.2f" % fit["chi2r_sobrecosto"]) if "chi2r_sobrecosto" in fit else "-",
            ("%.1f" % fit["p_max"]) if "p_max" in fit else "-", ("%.2f" % fit["S_max"]) if "S_max" in fit else "-"))
    print()
    print("Fraccion serial de Karp-Flatt e(p) (constante = Amdahl puro; creciente = sobrecosto que crece con p)")
    print("%-22s " % "serie" + " ".join("%7s" % ("p=%d" % p) for p in P_LISTA[1:]))
    for nombre, fit in series:
        if not fit:
            continue
        kf = dict((p, e) for p, e, se in fit["karp_flatt"])
        print("%-22s " % nombre + " ".join("%7s" % (("%.3f" % kf[p]) if p in kf else "-") for p in P_LISTA[1:]))


def git(*args):
    try:
        return subprocess.check_output(("git",) + args, cwd=REPO, stderr=subprocess.STDOUT).decode().strip()
    except (OSError, subprocess.CalledProcessError) as e:
        out = getattr(e, "output", b"")
        return "ERROR: " + (out.decode().strip() if out else str(e))


def main():
    ap = argparse.ArgumentParser(description="Graficas y ajustes de la Tarea A.3")
    ap.add_argument("--guardar", action="store_true", help="hace commit de las figuras y las tablas")
    a = ap.parse_args()
    a1 = leer_csv("parteA1_resumen.csv")
    a2 = leer_csv("parteA2_resumen.csv")
    comp = leer_csv("parteA2_mpi_vs_omp.csv")
    faltan = [n for n, v in (("datos/parteA1_resumen.csv", a1), ("datos/parteA2_resumen.csv", a2),
                             ("datos/parteA2_mpi_vs_omp.csv", comp)) if v is None]
    if faltan:
        print("[FALLA] faltan " + ", ".join(faltan) + ": corre antes  python3 analisis/parteA1_resumen.py  y  "
              "python3 analisis/parteA2_resumen.py")
        sys.exit(1)
    mpi = serie(a1, "S", "S_desv")
    omp = serie(a2, "S", "S_desv")
    mpi_i = serie(comp, "S_iter_mpi", "S_iter_mpi_desv")
    omp_i = serie(comp, "S_iter_omp", "S_iter_omp_desv")
    fmpi, fomp, fmpi_i, fomp_i = ajustar(mpi), ajustar(omp), ajustar(mpi_i), ajustar(omp_i)
    series = [("MPI corrida", fmpi), ("OpenMP corrida", fomp), ("MPI por iteracion", fmpi_i), ("OpenMP por iteracion", fomp_i)]
    imprimir(series)
    tablas(series)
    print()
    print("Tablas guardadas en datos/parteA3_ajustes.csv, parteA3_ajustes.md y parteA3_karp_flatt.csv")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print()
        print("[AVISO] falta matplotlib, asi que no se hicieron las figuras. Instalalo una vez con")
        print("    python3 -m pip install --user matplotlib")
        print("y repite este comando.")
        sys.exit(1)
    estilo(plt)
    hechas = figuras(plt, mpi, omp, fmpi, fomp, mpi_i, omp_i, fmpi_i, fomp_i)
    print("Figuras guardadas: " + ", ".join(os.path.relpath(h, REPO) for h in hechas))
    if a.guardar:
        out = git("add", "datos", "reporte/figuras", "analisis/parteA3_graficas.py")
        if out.startswith("ERROR"):
            print("AVISO: git add fallo: " + out)
            return
        if git("diff", "--cached", "--quiet").startswith("ERROR"):
            out = git("commit", "-q", "-m", "Tarea A.3: graficas de S(p) y E(p), ajustes de Amdahl y fraccion de Karp-Flatt")
            print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
        else:
            print("Nada nuevo que guardar en Git.")


if __name__ == "__main__":
    main()
