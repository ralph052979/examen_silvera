#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reporte_tablas.py - escribe las tablas del reporte (reporte/tablas/*.tex) a partir de las
tablas de datos/ que generan los demas scripts de analisis/. Asi ningun numero de las tablas
del reporte se copia a mano.

Uso (desde la carpeta del repositorio, despues de los scripts de las Partes A-D):
  python3 analisis/reporte_tablas.py
Solo usa la biblioteca estandar. Todo lo que imprime es ASCII.
"""
from __future__ import print_function, division
import csv
import io
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(REPO, "datos")
SALIDA = os.path.join(REPO, "reporte", "tablas")


def leer(nombre):
    ruta = os.path.join(DATOS, nombre)
    if not os.path.isfile(ruta):
        print("[FALLA] falta datos/%s: corre antes el script de analisis que la genera" % nombre)
        sys.exit(1)
    with io.open(ruta, "r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def x(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def pm(m, s, f="%.1f"):
    if m is None:
        return "--"
    return (f + r" $\pm$ " + f) % (m, s or 0.0)


def nf(v, f="%.2f", div=1.0):
    """numero con formato, o -- si falta"""
    v = x(v) if not isinstance(v, float) else v
    return "--" if v is None else f % (v / div)


def escribir(nombre, lineas):
    os.path.isdir(SALIDA) or os.makedirs(SALIDA)
    with io.open(os.path.join(SALIDA, nombre), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(lineas) + u"\n")
    print("  reporte/tablas/" + nombre)


def sci(v, cifras=1):
    """1.2e-08 -> $1.2\\times10^{-8}$"""
    if v is None:
        return "--"
    if v == 0:
        return "0"
    m, e = ("%.*e" % (cifras, v)).split("e")
    return r"$%s\times10^{%d}$" % (m, int(e))


def main():
    print("Tablas del reporte:")
    # ---------------- Tarea 0: caso base
    k7 = ([f for f in leer("tarea0_calibracion.csv") if f["cand"] == "K7"]
          if os.path.isfile(os.path.join(DATOS, "tarea0_calibracion.csv")) else [])
    a1p = [f for f in leer("parteA1_resumen.csv") if f["p"] == "1"]
    if k7:
        k = k7[0]
        t1 = (u"%.1f min (%.0f $\\pm$ %.0f s, A.1)" % (x(a1p[0]["t_media_s"]) / 60.0, x(a1p[0]["t_media_s"]),
                                                    x(a1p[0]["t_desv_s"]))) if a1p else u"%.1f min" % x(k["t1_min"])
        izq = [(u"Sistema", u"2H-WSe$_2$ en volumen, $2\\times2\\times1$ con V$_\\mathrm{Se}$"),
               (u"Átomos (W + Se)", u"%d (8 + 15)" % int(x(k["nions"]))),
               (u"Electrones de valencia", u"%d" % int(x(k["nelect"]))),
               (u"Funcional y PAW", u"PBE; W\\_sv, Se"),
               (u"ENCUT, PREC", u"%d eV, Normal" % int(x(k["encut"]))),
               (u"k-points", u"Gamma %s $\\to$ %d irreducibles" % (k["malla"].replace("x", u"$\\times$"), int(x(k["nkpts"]))))]
        der = [(u"NBANDS", u"%d" % int(x(k["nbands"]))),
               (u"Ondas planas por k-point (máx.)", u"%d" % int(x(k["pw_max"]))),
               (u"Malla FFT", k["fft"].replace(" x ", u"$\\times$")),
               (u"Iteraciones SCF (EDIFF = 1E-6 eV)", u"%d" % int(x(k["n_iter"]))),
               (u"Memoria estimada (VASP) / medida", u"%.0f MB / %.0f MB" % (x(k["mem_est_mb"]), x(k["mem_max_mb"]))),
               (u"$t$(1 núcleo)", t1)]
        L = [r"\begin{tabular}{ll@{\hspace{1.5em}}ll}", r"\toprule"]
        for (a1_, b1_), (a2_, b2_) in zip(izq, der):
            L.append(u"%s & %s & %s & %s \\\\" % (a1_, b1_, a2_, b2_))
        L += [r"\bottomrule", r"\end{tabular}"]
        escribir("tab_T0.tex", L)
    # ---------------- A.1 y A.2
    a1 = dict((int(x(f["p"])), f) for f in leer("parteA1_resumen.csv"))
    a2 = dict((int(x(f["p"])), f) for f in leer("parteA2_resumen.csv"))
    it = dict((int(x(f["p"])), f) for f in leer("parteA2_mpi_vs_omp.csv"))
    L = [r"\begin{tabular}{r rrr rrr rr}", r"\toprule",
         r" & \multicolumn{3}{c}{MPI (A.1): SCF completo} & \multicolumn{3}{c}{OpenMP (A.2): 4 iteraciones} & \multicolumn{2}{c}{s por iteración SCF} \\",
         r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-9}",
         r"$p$ & $t$ (s) & $S$ & $E$ & $t$ (s) & $S$ & $E$ & MPI & OpenMP \\", r"\midrule"]
    for p in sorted(a1):
        f, g, h = a1[p], a2.get(p), it.get(p)
        L.append(r"%d & %s & %.2f & %.3f & %s & %.2f & %.3f & %.1f & %.1f \\" % (
            p, pm(x(f["t_media_s"]), x(f["t_desv_s"])), x(f["S"]), x(f["E"]),
            pm(x(g["t_media_s"]), x(g["t_desv_s"])), x(g["S"]), x(g["E"]),
            x(h["t_iter_mpi_s"]), x(h["t_iter_omp_s"])))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_A.tex", L)
    # ---------------- A.3 ajustes
    aj = leer("parteA3_ajustes.csv")
    nombres = {"MPI corrida": "MPI, corrida", "OpenMP corrida": "OpenMP, corrida",
               "MPI por iteracion": "MPI, por iteración", "OpenMP por iteracion": "OpenMP, por iteración"}
    kf = {}
    for f in leer("parteA3_karp_flatt.csv"):
        kf.setdefault(f["serie"], {})[int(x(f["p"]))] = x(f["e"])
    L = [r"\begin{tabular}{l rr r rrr rr}", r"\toprule",
         r" & \multicolumn{3}{c}{Amdahl} & \multicolumn{3}{c}{con sobrecosto $c\,(p-1)$} & \multicolumn{2}{c}{Karp--Flatt $e(p)$} \\",
         r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-9}",
         r"serie & $f$ & $\chi^2_\nu$ & $1/f$ & $f$ & $c$ & $p^{*}$ & $p=2$ & $p=24$ \\", r"\midrule"]
    for f in aj:
        k = kf.get(f["serie"], {})
        L.append(r"%s & %s & %.0f & %.1f & %s & %s & %s & %.3f & %.3f \\" % (
            nombres.get(f["serie"], f["serie"]), pm(x(f["f"]), x(f["f_desv"]), "%.3f"), x(f["chi2r_amdahl"]),
            x(f["S_inf_amdahl"]), pm(x(f["f2"]), x(f["f2_desv"]), "%.3f"), sci(x(f["c"])),
            ("%.1f" % x(f["p_max"])) if x(f["p_max"]) else "--", k.get(2, float("nan")), k.get(24, float("nan"))))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_A3.tex", L)
    # ---------------- B.1
    b1 = leer("parteB1_resumen.csv")
    rut = dict((int(x(f["p"])), f) for f in leer("parteB1_rutinas.csv"))
    L = [r"\begin{tabular}{r l r r r r r r r}", r"\toprule",
         r"$p$ & malla & NKPTS & $t$ (s) & $E_w$ & s/iter. & EDDAV (s/iter.) & mem./proceso (GB) & mem. total (GB) \\",
         r"\midrule"]
    for f in b1:
        p = int(x(f["p"]))
        ed = x(rut[p]["EDDAV_s"]) if p in rut else None
        L.append(r"%d & %s & %d & %s & %s & %s & %s & %s & %s \\" % (
            p, f["malla"].replace("x", r"$\times$"), int(x(f["nkpts"])), pm(x(f["t_media_s"]), x(f["t_desv_s"])),
            pm(x(f["Ew"]), x(f["Ew_desv"]), "%.3f"), nf(f["t_iter_s"]), "--" if ed is None else "%.2f" % ed,
            nf(f["mem_proc_mb"], div=1024.0), nf(f["mem_total_mb"], div=1024.0)))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_B1.tex", L)
    # ---------------- B.2
    b2 = leer("parteB2_resumen.csv")
    L = [r"\begin{tabular}{c r r r r r r c}", r"\toprule",
         r"$N_\mathrm{MPI}\times N_\mathrm{OMP}$ & $n$ & $t$ (s) & s/iter. & $S=t(1)/t$ & mem./proceso (GB) & mem. total (GB) & 2 sockets \\",
         r"\midrule"]
    for f in b2:
        dos = [v for v in (f.get("en_2_sockets") or "").split("/") if v not in ("", "-")]
        L.append(r"%d$\times$%d & %d & %s & %s & %s & %s & %s & %s \\" % (
            int(x(f["n_mpi"])), int(x(f["n_omp"])), int(x(f["n"])), pm(x(f["t_media_s"]), x(f["t_desv_s"])),
            nf(f["t_iter_s"], "%.1f"), pm(x(f["S"]), x(f["S_desv"]), "%.2f"), nf(f["mem_proc_mb"], div=1024.0),
            nf(f["mem_total_mb"], div=1024.0), ("%s de %d" % (dos[0], int(x(f["n_mpi"])))) if dos else "--"))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_B2.tex", L)
    # ---------------- C.1
    c1 = leer("parteC1_resumen.csv")
    L = [r"\begin{tabular}{r c l r r r}", r"\toprule",
         r"$p$ & valores distintos en 3 rep. & $E(p)$ (eV, OSZICAR) & $E(p)-E(1)$ (eV) & $|\Delta E|/|E(1)|$ & $|\Delta E|/(\varepsilon_\mathrm{mach}|E(1)|)$ \\",
         r"\midrule"]
    for f in c1:
        e = float(f["e_oszicar_eV"])
        L.append(r"%d & %s & %.9f & %s & %s & %s \\" % (
            int(x(f["p"])), f["valores_distintos"], e, sci(x(f["de_eV"])), sci(x(f["de_relativa"])), sci(x(f["n_eps"]))))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_C1.tex", L)
    it4 = leer("parteC1_iteracion4.csv")
    grupos = {}
    for f in it4:
        g = grupos.setdefault(int(x(f["n_mpi"])), {"T": [], "corridas": 0, "vals": set()})
        g["T"].append(int(x(f["n_omp"])))
        g["corridas"] += int(x(f["corridas"]))
        g["vals"].update(float(v) for v in f["e_iter4_eV"].split(" | "))
    ref = sorted(grupos[1]["vals"])[0] if 1 in grupos else None
    L = [r"\begin{tabular}{r l r l r}", r"\toprule",
         r"$N_\mathrm{MPI}$ & $N_\mathrm{OMP}$ medidos & corridas & $E$ tras 4 iteraciones (eV) & $E - E(1\times1)$ (eV) \\",
         r"\midrule"]
    for n in sorted(grupos):
        g = grupos[n]
        vals = sorted(g["vals"], key=lambda v: abs(v - ref) if ref is not None else 0.0)
        e = u" / ".join(u"%.9f" % v for v in vals)
        de = vals[0] - ref if ref is not None else None
        L.append(r"%d & %s & %d & %s & %s \\" % (n, ", ".join(str(t) for t in sorted(g["T"])), g["corridas"], e,
                                               "0" if not de else sci(de)))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_C1_iter4.tex", L)
    red = leer("parteC1_reduccion.csv")
    L = [r"\begin{tabular}{l r r r}", r"\toprule",
         r"precisión & $p$ & diferencia con $p=1$ ($\varepsilon_\mathrm{mach}$) & error frente a la suma exacta ($\varepsilon_\mathrm{mach}$) \\",
         r"\midrule"]
    for f in red:
        L.append(r"%s & %d & %.0f & %.0f \\" % (f["precision"].replace("(", "(\\texttt{").replace(")", "})"), int(x(f["p"])),
                                                x(f["dif_vs_p1_eps"]), x(f["error_vs_exacta_eps"])))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_C1_reduccion.tex", L)
    # ---------------- C.2
    c2 = leer("parteC2_resumen.csv")
    L = [r"\begin{tabular}{r r r r l r r}", r"\toprule",
         r"ENCUT (eV) & $t$ (s) & $t/t(400)$ & ondas planas & malla FFT & mem. (MB) & error (meV/átomo) \\", r"\midrule"]
    for f in c2:
        err = x(f["error_meV_atomo"])
        L.append(r"%d & %s & %s & %s & %s & %s & %s \\" % (
            int(x(f["encut_eV"])), pm(x(f["t_media_s"]), x(f["t_desv_s"])), nf(f["t_rel_400"]), nf(f["npw"], "%.0f"),
            f["fft"].replace("x", r"$\times$") if f.get("fft") else "--", nf(f["mem_mb"], "%.0f"),
            "referencia" if err is None else "%.3f" % err))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_C2.tex", L)
    # ---------------- D
    co = leer("parteD_costo.csv")
    L = [r"\begin{tabular}{l r r r}", r"\toprule",
         r"configuración (SCF completo del caso base) & núcleos & $t$ (min) & h-núcleo por cálculo \\", r"\midrule"]
    for f in co:
        conf = f["configuracion"].replace(" (A.1)", "").replace(" (B.1, malla del caso base)", " (B.1)")
        conf = conf.replace("1 procesos MPI", "1 proceso MPI") if conf.startswith("1 procesos") else conf
        L.append(r"%s & %d & %.1f & %.2f \\" % (conf, int(x(f["nucleos"])), x(f["t_s"]) / 60.0, x(f["h_nucleo"])))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_D_costo.tex", L)
    cons = leer("parteD_consumo.csv")
    bonitos = {"Tarea 0: calibracion final (K1-K7)": u"Tarea 0: calibración final (K1--K7)",
               "A.2: diagnostico de OpenMP (2 rondas)": u"A.2: diagnóstico de OpenMP (2 rondas)",
               "A.2: barrido OpenMP (4 iteraciones)": u"A.2: barrido OpenMP (4 iteraciones)",
               "B.1: escalamiento debil": u"B.1: escalamiento débil",
               "B.2: hibrido MPI x OpenMP": u"B.2: híbrido MPI $\\times$ OpenMP",
               "Pruebas e intentos previos (17-18 nov)": u"Pruebas e intentos previos (17--18 nov.)"}
    L = [r"\begin{tabular}{l r r}", r"\toprule", r"parte & jobs & h-núcleo \\", r"\midrule"]
    for f in cons:
        if f["categoria"] == "TOTAL":
            L.append(r"\midrule")
            L.append(r"\textbf{total (de 250)} & %s & \textbf{%.2f} \\" % (f["jobs"], x(f["h_nucleo"])))
        else:
            L.append(r"%s & %s & %.2f \\" % (bonitos.get(f["categoria"], f["categoria"]), f["jobs"], x(f["h_nucleo"])))
    L += [r"\bottomrule", r"\end{tabular}"]
    escribir("tab_D_consumo.tex", L)


if __name__ == "__main__":
    main()
