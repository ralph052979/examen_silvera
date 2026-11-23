#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteB_graficas.py - figuras de la Parte B (escalamiento debil e hibrido).

Lee las tablas de analisis/parteB_resumen.py (datos/parteB1_resumen.csv y parteB2_resumen.csv)
y los OUTCAR de las corridas de la B.1 (logs/jobs_B1.txt y corridas/) para separar el tiempo
de cada iteracion SCF por rutina de VASP.

Escribe:
  reporte/figuras/B1_debil.pdf y .png     (a) Ew(p) de la corrida, por iteracion y del
                                           eigensolver (EDDAV); (b) memoria por proceso y
                                           total frente a p, con el limite de la QOS
  reporte/figuras/B2_hibrido.pdf y .png   (a) tiempo por combinacion N_MPI x N_OMP;
                                           (b) memoria total y por proceso
  datos/parteB1_rutinas.csv               tiempo medio por iteracion de cada rutina (B.1)
  datos/parteB1_memoria.md                modelo de memoria: m(p) = a + b NKPTS y el p
                                           a partir del cual el total pasa el limite de la QOS

Uso (desde la carpeta del repositorio):
  python3 analisis/parteB_graficas.py             # figuras y tablas
  python3 analisis/parteB_graficas.py --guardar   # ademas hace commit
Necesita matplotlib. Todo lo que imprime es ASCII.
"""
from __future__ import print_function, division
import argparse
import csv
import glob
import io
import math
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(REPO, "datos")
FIGS = os.path.join(REPO, "reporte", "figuras")
P_LISTA = [1, 2, 4, 8, 12, 16, 24]
QOS_MB = 62600.0
RUTINAS = ["EDDAV", "POTLOK", "CHARGE", "SETDIJ", "MIXING", "DOS", "LOOP"]
C_AZUL, C_NARANJA, C_AQUA = "#2a78d6", "#eb6834", "#1baf7a"
TINTA, TINTA2, MUTED, REJILLA, EJE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"


def leer(ruta):
    try:
        with io.open(ruta, "r", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except (IOError, OSError):
        return ""


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


def git(*args):
    try:
        return subprocess.check_output(("git",) + args, cwd=REPO, stderr=subprocess.STDOUT).decode().strip()
    except (OSError, subprocess.CalledProcessError) as e:
        out = getattr(e, "output", b"")
        return "ERROR: " + (out.decode().strip() if out else str(e))


def rutinas_b1():
    """tiempo medio por iteracion de cada rutina, por p (media de las repeticiones)"""
    acum = {}
    for l in leer(os.path.join(REPO, "logs", "jobs_B1.txt")).splitlines():
        c = l.split()
        if len(c) < 3 or not c[2].isdigit():
            continue
        P, R, J = int(c[0]), c[1], c[2]
        d = os.path.join(REPO, "corridas", "ex_B1_p%d_r%s_%s" % (P, R, J))
        t = leer(os.path.join(d, "OUTCAR"))
        if not t:
            continue
        for r in RUTINAS:
            vals = [float(v) for v in re.findall(r"^\s+%s:\s+cpu time\s+[\d.]+:\s+real time\s+([\d.]+)" % r, t, re.M)]
            if vals:
                acum.setdefault((P, r), []).append(sum(vals) / len(vals))
    filas = []
    for P in P_LISTA:
        f = {"p": P}
        for r in RUTINAS:
            v = acum.get((P, r))
            f[r] = sum(v) / len(v) if v else None
        if any(f[r] is not None for r in RUTINAS):
            filas.append(f)
    return filas


def estilo(plt):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.edgecolor": EJE, "axes.linewidth": 0.8, "axes.labelcolor": TINTA, "xtick.color": TINTA2,
        "ytick.color": TINTA2, "text.color": TINTA, "grid.color": REJILLA, "grid.linewidth": 0.6, "axes.grid": True,
        "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "pdf.fonttype": 42, "savefig.dpi": 200, "hatch.color": "white", "hatch.linewidth": 1.0})


def eje_p(ax, hasta=30):
    try:
        ax.set_xscale("log", base=2)
    except (TypeError, ValueError):
        ax.set_xscale("log", basex=2)
    ax.set_xticks(P_LISTA)
    ax.set_xticklabels([str(p) for p in P_LISTA])
    ax.minorticks_off()
    ax.set_xlim(0.85, hasta)
    ax.set_xlabel(u"p (procesos MPI = núcleos)")


def guardar(fig, nombre):
    hechas = []
    os.path.isdir(FIGS) or os.makedirs(FIGS)
    for ext in ("pdf", "png"):
        ruta = os.path.join(FIGS, nombre + "." + ext)
        fig.savefig(ruta, bbox_inches="tight")
        hechas.append(ruta)
    return hechas


def modelo_memoria(b1):
    """m(p) = a + b NKPTS por proceso (minimos cuadrados) y p con p m(p) = limite de la QOS"""
    pts = [(x(f["nkpts"]), x(f["mem_proc_mb"])) for f in b1 if x(f.get("mem_proc_mb"))]
    n = len(pts)
    if n < 2:
        return None
    sx = sum(a for a, b in pts); sy = sum(b for a, b in pts)
    sxx = sum(a * a for a, b in pts); sxy = sum(a * b for a, b in pts)
    pend = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    ordn = (sy - pend * sx) / n
    # total(p) = p (a + 3 b p) = QOS  ->  3 b p^2 + a p - QOS = 0
    A, B, C = 3.0 * pend, ordn, -QOS_MB
    disc = B * B - 4 * A * C
    p_lim = (-B + math.sqrt(disc)) / (2 * A) if (A > 0 and disc >= 0) else None
    res = max(abs(b - (ordn + pend * a)) for a, b in pts)
    return {"a_mb": ordn, "b_mb_por_kpt": pend, "p_limite": p_lim, "residuo_max_mb": res}


def fig_b1(plt, b1, rut, mod):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.4))
    ps = [int(x(f["p"])) for f in b1]
    a1.plot([1, 24], [1, 1], "--", color=MUTED, lw=1.0, label="ideal", zorder=2)
    a1.errorbar(ps, [x(f["Ew"]) for f in b1], yerr=[x(f["Ew_desv"]) or 0 for f in b1], fmt="o-", color=C_AZUL, ms=5.5,
                lw=1.6, mec="white", mew=0.8, capsize=2, elinewidth=0.9, label=u"corrida completa", zorder=4)
    a1.plot(ps, [x(f["Ew_iter"]) for f in b1], "s--", color=C_NARANJA, ms=5, lw=1.4, mec="white", mew=0.8,
            label=u"por iteración SCF", zorder=3)
    if rut and all(f["EDDAV"] for f in rut):
        e1 = rut[0]["EDDAV"]
        a1.plot([f["p"] for f in rut], [e1 / f["EDDAV"] for f in rut], "^:", color=C_AQUA, ms=5, lw=1.6, mec="white",
                mew=0.8, label=u"solo el eigensolver (EDDAV)", zorder=3)
    eje_p(a1)
    a1.set_ylim(0, 1.1)
    a1.set_ylabel(u"eficiencia débil E$_w$(p) = t(1)/t(p)")
    a1.set_title("(a)", loc="left", fontsize=9, fontweight="bold")
    a1.legend(loc="lower left")
    gb = 1024.0
    a2.plot(ps, [x(f["mem_total_mb"]) / gb for f in b1], "o-", color=C_AZUL, ms=5.5, lw=1.6, mec="white", mew=0.8,
            label=u"total (suma)", zorder=4)
    a2.plot(ps, [x(f["mem_proc_mb"]) / gb for f in b1], "s-", color=C_NARANJA, ms=5, lw=1.6, mec="white", mew=0.8,
            label=u"por proceso (máx.)", zorder=4)
    if mod:
        pp = [1 + 0.25 * i for i in range(int((32 - 1) / 0.25) + 1)]
        a2.plot(pp, [p * (mod["a_mb"] + 3 * mod["b_mb_por_kpt"] * p) / gb for p in pp], ":", color=C_AZUL, lw=1.3,
                zorder=3, label=u"modelo p·(a + b·3p)")
    a2.axhline(QOS_MB / gb, ls="--", lw=1.0, color=MUTED, zorder=2)
    a2.text(1.0, QOS_MB / gb * 1.12, u"límite de memoria de la QOS (62 600 MB)", fontsize=7.5, color=TINTA2, va="bottom")
    if mod and mod["p_limite"]:
        a2.axvline(mod["p_limite"], ls=":", lw=0.9, color=MUTED, zorder=1)
        a2.text(mod["p_limite"] * 0.96, 2.0, u"p ≈ %.0f" % mod["p_limite"], fontsize=7.5, color=TINTA2, ha="right")
    eje_p(a2, hasta=34)
    a2.set_yscale("log")
    a2.set_yticks([0.1, 0.3, 1, 3, 10, 30, 100])
    a2.set_yticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])
    a2.minorticks_off()
    a2.set_ylim(0.1, 150)
    a2.set_ylabel(u"memoria máxima (GB)")
    a2.set_title("(b)", loc="left", fontsize=9, fontweight="bold")
    a2.legend(loc="upper left", bbox_to_anchor=(0.0, 0.88))
    fig.tight_layout()
    h = guardar(fig, "B1_debil")
    plt.close(fig)
    return h


def fig_b2(plt, b2):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3))
    et = [u"%d×%d" % (int(x(f["n_mpi"])), int(x(f["n_omp"]))) for f in b2]
    xs = list(range(len(b2)))
    dos = [any(v not in ("0", "-", "") for v in (f.get("en_2_sockets") or "-").split("/")) for f in b2]
    t = [x(f["t_media_s"]) for f in b2]
    barras = a1.bar(xs, t, yerr=[x(f["t_desv_s"]) or 0 for f in b2], width=0.72, color=C_AZUL, edgecolor="white",
                    linewidth=0.8, capsize=2, error_kw={"elinewidth": 0.9, "ecolor": TINTA2}, zorder=3)
    for b, d in zip(barras, dos):
        if d:
            b.set_hatch("////")
    imin = min(range(len(t)), key=lambda i: t[i])
    a1.annotate(u"óptimo\n%.0f s" % t[imin], (xs[imin], t[imin]), xytext=(0, 14), textcoords="offset points",
                ha="center", fontsize=7.5, color=TINTA2)
    a1.annotate(u"%.0f s" % t[-1], (xs[-1], t[-1]), xytext=(0, 4), textcoords="offset points", ha="center",
                fontsize=7.5, color=TINTA2)
    a1.set_xticks(xs)
    a1.set_xticklabels(et)
    a1.set_xlabel(u"N$_{MPI}$ × N$_{OMP}$ (24 núcleos)")
    a1.set_ylabel(u"tiempo de pared, 4 iteraciones SCF (s)")
    a1.set_title("(a)", loc="left", fontsize=9, fontweight="bold")
    a1.grid(axis="x", visible=False)
    from matplotlib.patches import Patch
    a1.legend(handles=[Patch(facecolor=C_AZUL, edgecolor="white", hatch="////",
                             label=u"algún proceso abarca los 2 sockets")], loc="upper left")
    w = 0.38
    a2.bar([i - w / 2 for i in xs], [x(f["mem_total_mb"]) / 1024.0 for f in b2], width=w, color=C_AZUL,
           edgecolor="white", linewidth=0.8, label=u"total", zorder=3)
    a2.bar([i + w / 2 for i in xs], [x(f["mem_proc_mb"]) / 1024.0 for f in b2], width=w, color=C_NARANJA,
           edgecolor="white", linewidth=0.8, label=u"por proceso (máximo)", zorder=3)
    a2.set_xticks(xs)
    a2.set_xticklabels(et)
    a2.set_xlabel(u"N$_{MPI}$ × N$_{OMP}$ (24 núcleos)")
    a2.set_ylabel(u"memoria máxima (GB)")
    a2.set_title("(b)", loc="left", fontsize=9, fontweight="bold")
    a2.grid(axis="x", visible=False)
    a2.legend(loc="upper right")
    fig.tight_layout()
    h = guardar(fig, "B2_hibrido")
    plt.close(fig)
    return h


def main():
    ap = argparse.ArgumentParser(description="Figuras de la Parte B")
    ap.add_argument("--guardar", action="store_true", help="hace commit de las figuras y tablas")
    a = ap.parse_args()
    b1 = leer_csv("parteB1_resumen.csv")
    b2 = leer_csv("parteB2_resumen.csv")
    if not b1 or not b2:
        print("[FALLA] faltan datos/parteB1_resumen.csv o parteB2_resumen.csv: corre antes  python3 analisis/parteB_resumen.py")
        sys.exit(1)
    rut = rutinas_b1()
    mod = modelo_memoria(b1)
    print("B.1: tiempo medio por iteracion SCF de cada rutina de VASP (s), media de las repeticiones")
    print("%-4s " % "p" + " ".join("%8s" % r for r in RUTINAS))
    for f in rut:
        print("%-4d " % f["p"] + " ".join("%8s" % ("-" if f[r] is None else "%.2f" % f[r]) for r in RUTINAS))
    with io.open(os.path.join(DATOS, "parteB1_rutinas.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["p"] + ["%s_s" % r for r in RUTINAS])
        for f in rut:
            w.writerow([f["p"]] + ["" if f[r] is None else "%.4f" % f[r] for r in RUTINAS])
    if mod:
        print()
        print("Memoria por proceso en la B.1: m = %.1f MB + %.2f MB x NKPTS (residuo maximo %.0f MB)" % (
            mod["a_mb"], mod["b_mb_por_kpt"], mod["residuo_max_mb"]))
        if mod["p_limite"]:
            print("Con NKPTS = 3p, la memoria total p m(p) pasa el limite de la QOS (62 600 MB) en p = %.1f" % mod["p_limite"])
        else:
            print("El modelo lineal de memoria no cruza el limite de la QOS (pendiente no positiva)")
        md = [u"# B.1 - modelo de memoria", u"",
              u"Memoria maxima por proceso (VmHWM) frente al numero total de k-points: m = a + b NKPTS, por minimos cuadrados "
              u"con las 7 medias de la B.1.", u"",
              u"| a (MB) | b (MB por k-point) | residuo maximo (MB) | p con p m(p) = 62 600 MB (NKPTS = 3p) |",
              u"|---|---|---|---|",
              u"| %.1f | %.2f | %.0f | %s |" % (mod["a_mb"], mod["b_mb_por_kpt"], mod["residuo_max_mb"],
                                              ("%.1f" % mod["p_limite"]) if mod["p_limite"] else u"-"), u""]
        with io.open(os.path.join(DATOS, "parteB1_memoria.md"), "w", encoding="utf-8") as fh:
            fh.write(u"\n".join(md))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[AVISO] falta matplotlib (python3 -m pip install --user matplotlib)")
        sys.exit(1)
    estilo(plt)
    hechas = fig_b1(plt, b1, rut, mod) + fig_b2(plt, b2)
    print()
    print("Figuras guardadas: " + ", ".join(os.path.relpath(h, REPO) for h in hechas))
    if a.guardar:
        out = git("add", "datos", "reporte/figuras", "analisis/parteB_graficas.py")
        if out.startswith("ERROR"):
            print("AVISO: git add fallo: " + out)
            return
        if git("diff", "--cached", "--quiet").startswith("ERROR"):
            out = git("commit", "-q", "-m", "Parte B: figuras del escalamiento debil y del hibrido, rutinas y modelo de memoria")
            print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
        else:
            print("Nada nuevo que guardar en Git.")


if __name__ == "__main__":
    main()
