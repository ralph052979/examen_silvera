#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteC_analisis.py - Parte C: reproducibilidad numerica (C.1) y precision frente a costo (C.2).

C.1  Energia total del caso base con p = 1, 2, 4, 8, 12, 16, 24 procesos MPI (corridas de la A.1:
     SCF completo, 3 repeticiones por p). Para cada p: la energia final, si las 3 repeticiones dan
     el mismo numero, y la diferencia respecto de p = 1 en eV y en unidades de eps_mach:
         n_eps = |E(p) - E(1)| / (eps_mach |E(1)|),   eps_mach = 2^-52 (doble precision).
     Fuente principal: la ultima linea DAV del OSZICAR (12 cifras significativas: resolucion de
     1E-9 eV); control: energy(sigma->0) del OUTCAR (8 decimales).
     Ademas compara la energia despues de exactamente 4 iteraciones (sin criterio de parada) con
     N procesos MPI (A.1 y B.2) y con T hilos OpenMP (A.2 y B.2): separa el efecto del numero de
     procesos MPI (reparto de las bandas) del de los hilos.
     Demostracion del orden de las reducciones: suma de 1 048 576 numeros aleatorios en p bloques
     (como una reduccion con p procesos), en precision simple (float32) y doble (float64).
C.2  Serie de ENCUT (malla 3x3x1, 1 nucleo, SCF completo, 3 repeticiones): error por atomo
     |E(ENCUT) - E(800 eV)| / N_atomos frente al tiempo de pared, y crecimiento relativo del tiempo,
     del numero de ondas planas y de la malla FFT. Tambien la convergencia con la malla de k-points
     (Tarea 0 y B.1, sin corridas nuevas).

Lee:     logs/jobs_A1.txt, logs/jobs_A2.txt, logs/jobs_B2.txt, logs/jobs_C2.txt, logs/jobs_t0.txt,
         logs/jobs_B1.txt y las carpetas corridas/<job>/ (OSZICAR, OUTCAR, tiempo.txt)
Escribe: datos/parteC1_energias.csv, parteC1_resumen.csv, parteC1_iteracion4.csv,
         parteC1_reduccion.csv, parteC2_corridas.csv, parteC2_resumen.csv, parteC2_kpuntos.csv,
         parteC_tabla.md y reporte/figuras/C2_precision_costo.pdf y .png

Uso (desde la carpeta del repositorio):
  python3 analisis/parteC_analisis.py             # tablas y figura
  python3 analisis/parteC_analisis.py --guardar   # ademas hace commit de tablas y figura

Necesita numpy para la demostracion de la reduccion y matplotlib para la figura (las tablas de
energia y tiempo se escriben igual). Todo lo que imprime es ASCII.
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
EPS64 = 2.0 ** -52
EPS32 = 2.0 ** -23
P_C1 = [1, 2, 4, 8, 12, 16, 24]
ENCUTS = [250, 300, 350, 400, 450, 500, 550, 600, 700, 800]
E_REF = 800
EDIFF = 1e-6
C_AZUL, C_NARANJA, C_AQUA = "#2a78d6", "#eb6834", "#1baf7a"
TINTA, TINTA2, MUTED, REJILLA, EJE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"


# ----------------------------------------------------------------------------- utilidades
def leer(ruta):
    try:
        with io.open(ruta, "r", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except (IOError, OSError):
        return ""


def clave_valor(ruta):
    d = {}
    for l in leer(ruta).splitlines():
        if "=" in l:
            k, v = l.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def media_desv(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return None, None
    m = sum(xs) / float(len(xs))
    if len(xs) == 1:
        return m, 0.0
    return m, math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def registro(nombre, ncampos):
    """lineas de logs/<nombre>: ncampos valores (texto) antes del jobid"""
    filas = []
    for l in leer(os.path.join(REPO, "logs", nombre)).splitlines():
        p = l.split()
        if len(p) >= ncampos + 1 and p[ncampos].isdigit():
            filas.append((p[:ncampos], p[ncampos], p[ncampos + 1:]))
    return filas


def carpeta(patron):
    d = sorted(glob.glob(os.path.join(REPO, "corridas", patron)))
    return d[-1] if d else None


def energias_oszicar(d):
    """lista de energias (texto, 12 cifras) de las lineas DAV/RMM del OSZICAR"""
    return re.findall(r"^(?:DAV|RMM):\s+\d+\s+(\S+)", leer(os.path.join(d, "OSZICAR")), re.M)


def datos_outcar(d):
    t = leer(os.path.join(d, "OUTCAR"))
    o = {"texto": t}
    m = re.findall(r"energy\(sigma->0\)\s*=\s*([-\d.]+)", t)
    o["e_sigma0"] = m[-1] if m else None
    m = re.search(r"NIONS\s*=\s*(\d+)", t)
    o["nions"] = int(m.group(1)) if m else None
    m = re.search(r"maximum number of plane-waves:\s*(\d+)", t)
    o["npw"] = int(m.group(1)) if m else None
    m = re.search(r"dimension x,y,z NGX =\s*(\d+)\s+NGY =\s*(\d+)\s+NGZ =\s*(\d+)", t)
    o["fft"] = tuple(int(x) for x in m.groups()) if m else None
    m = re.search(r"Maximum memory used \(kb\):\s*([\d.]+)", t)
    o["mem_mb"] = float(m.group(1)) / 1024.0 if m else None
    o["convergido"] = "aborting loop because EDIFF is reached" in t
    o["n_iter"] = len(re.findall(r"LOOP:\s+cpu time", t))
    m = re.search(r"running\s+(\d+)\s+mpi-ranks,\s+with\s+(\d+)\s+threads/rank", t)
    o["ranks"], o["hilos"] = (int(m.group(1)), int(m.group(2))) if m else (None, None)
    m = re.search(r"NKPTS\s*=\s*(\d+)", t)
    o["nkpts"] = int(m.group(1)) if m else None
    return o


def git(*args):
    try:
        return subprocess.check_output(("git",) + args, cwd=REPO, stderr=subprocess.STDOUT).decode().strip()
    except (OSError, subprocess.CalledProcessError) as e:
        out = getattr(e, "output", b"")
        return "ERROR: " + (out.decode().strip() if out else str(e))


def escribir_csv(nombre, cols, filas):
    with io.open(os.path.join(DATOS, nombre), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for f in filas:
            w.writerow([("" if f.get(k) is None else (("%.12g" % f[k]) if isinstance(f[k], float) else f[k]))
                        for k in cols])


# ----------------------------------------------------------------------------- C.1
def c1_corridas():
    runs = []
    for (P, R), jid, _ in registro("jobs_A1.txt", 2):
        d = carpeta("ex_A1_p%s_r%s_%s" % (P, R, jid))
        if not d:
            continue
        es = energias_oszicar(d)
        o = datos_outcar(d)
        tie = clave_valor(os.path.join(d, "tiempo.txt"))
        runs.append({"serie": "A.1", "n_mpi": int(P), "n_omp": 1, "r": int(R), "jobid": jid,
                     "estado": tie.get("estado"), "n_iter": len(es), "convergido": int(o["convergido"]),
                     "e_final_oszicar": es[-1] if es else None, "e_iter4_oszicar": es[3] if len(es) >= 4 else None,
                     "e_sigma0_outcar": o["e_sigma0"], "ranks": o["ranks"], "hilos": o["hilos"]})
    for (P, R), jid, _ in registro("jobs_A2.txt", 2):
        d = carpeta("ex_A2_p%s_r%s_%s" % (P, R, jid))
        if not d:
            continue
        es = energias_oszicar(d)
        runs.append({"serie": "A.2", "n_mpi": 1, "n_omp": int(P), "r": int(R), "jobid": jid,
                     "estado": clave_valor(os.path.join(d, "tiempo.txt")).get("estado"), "n_iter": len(es),
                     "e_iter4_oszicar": es[3] if len(es) >= 4 else None})
    for (N, T, R), jid, _ in registro("jobs_B2.txt", 3):
        d = carpeta("ex_B2_n%sx%s_r%s_%s" % (N, T, R, jid))
        if not d:
            continue
        es = energias_oszicar(d)
        runs.append({"serie": "B.2", "n_mpi": int(N), "n_omp": int(T), "r": int(R), "jobid": jid,
                     "estado": clave_valor(os.path.join(d, "tiempo.txt")).get("estado"), "n_iter": len(es),
                     "e_iter4_oszicar": es[3] if len(es) >= 4 else None})
    return runs


def c1_resumen(runs):
    a1 = [c for c in runs if c["serie"] == "A.1" and c["estado"] == "OK" and c.get("e_final_oszicar")]
    base = [float(c["e_final_oszicar"]) for c in a1 if c["n_mpi"] == 1]
    base_out = [float(c["e_sigma0_outcar"]) for c in a1 if c["n_mpi"] == 1 and c.get("e_sigma0_outcar")]
    if not base:
        return [], None
    e1 = base[0]
    e1o = base_out[0] if base_out else None
    filas = []
    for P in P_C1:
        sel = [c for c in a1 if c["n_mpi"] == P]
        if not sel:
            continue
        vals = sorted(set(c["e_final_oszicar"] for c in sel))
        valso = sorted(set(c["e_sigma0_outcar"] for c in sel if c.get("e_sigma0_outcar")))
        e = float(sel[0]["e_final_oszicar"])
        de = e - e1
        deo = (float(valso[0]) - e1o) if (valso and e1o is not None) else None
        filas.append({"p": P, "n": len(sel), "iteraciones": "/".join(str(c["n_iter"]) for c in sel),
                      "valores_distintos": len(vals), "e_oszicar_eV": vals[0] if len(vals) == 1 else " | ".join(vals),
                      "de_eV": de, "de_relativa": abs(de) / abs(e1), "n_eps": abs(de) / (EPS64 * abs(e1)),
                      "e_outcar_eV": valso[0] if len(valso) == 1 else " | ".join(valso),
                      "de_outcar_eV": deo, "n_eps_outcar": (abs(deo) / (EPS64 * abs(e1o))) if deo is not None else None})
    return filas, e1


def c1_iteracion4(runs):
    """energia despues de 4 iteraciones por (N_MPI, N_OMP)"""
    grupos = {}
    for c in runs:
        if c.get("e_iter4_oszicar") and c.get("estado") == "OK":
            grupos.setdefault((c["n_mpi"], c["n_omp"]), []).append(c)
    ref = None
    for (n, t), sel in grupos.items():
        if (n, t) == (1, 1):
            ref = float(sel[0]["e_iter4_oszicar"])
    filas = []
    for (n, t) in sorted(grupos, key=lambda k: (k[1] != 1 and k[0] != 1, k[0] if k[1] == 1 else 100 + k[1], -k[0])):
        sel = grupos[(n, t)]
        vals = sorted(set(c["e_iter4_oszicar"] for c in sel))
        e = float(vals[0])
        filas.append({"n_mpi": n, "n_omp": t, "corridas": len(sel), "series": "+".join(sorted(set(c["serie"] for c in sel))),
                      "valores_distintos": len(vals), "e_iter4_eV": " | ".join(vals),
                      "de_eV": (e - ref) if ref is not None else None,
                      "n_eps": (abs(e - ref) / (EPS64 * abs(ref))) if ref is not None else None,
                      "bandas_por_proceso": 144 // n if 144 % n == 0 else None})
    return filas


def demo_reduccion():
    """suma de N numeros en p bloques contiguos (como una reduccion con p procesos)"""
    try:
        import numpy as np
    except ImportError:
        return None
    rng = np.random.RandomState(2026)
    n = 2 ** 20
    x64 = rng.uniform(0.0, 1.0, n)
    x32 = x64.astype(np.float32)
    exacta64 = math.fsum(x64.tolist())
    exacta32 = math.fsum(x32.astype(np.float64).tolist())
    filas = []
    for nombre, x, tipo, eps, exacta in (("doble (float64)", x64, np.float64, EPS64, exacta64),
                                         ("simple (float32)", x32, np.float32, EPS32, exacta32)):
        sumas = {}
        for p in (1, 4, 16, 24):
            bloques = np.array_split(x, p)
            parciales = [np.cumsum(b, dtype=tipo)[-1] for b in bloques]       # suma secuencial en cada proceso
            total = np.cumsum(np.array(parciales, dtype=tipo), dtype=tipo)[-1]  # reduccion de las parciales
            sumas[p] = float(total)
        for p in (1, 4, 16, 24):
            filas.append({"precision": nombre, "eps_mach": eps, "n_terminos": n, "p": p, "suma": sumas[p],
                          "dif_vs_p1_rel": abs(sumas[p] - sumas[1]) / abs(sumas[1]),
                          "dif_vs_p1_eps": abs(sumas[p] - sumas[1]) / (eps * abs(sumas[1])),
                          "error_vs_exacta_eps": abs(sumas[p] - exacta) / (eps * abs(exacta))})
    return filas


# ----------------------------------------------------------------------------- C.2
def c2_corridas():
    runs = []
    for (R,), jid, resto in registro("jobs_C2.txt", 1):
        lista = resto[0].split(":") if resto else []
        for E in lista:
            d = carpeta("ex_C2_r%s_%s_E%s" % (R, jid, E))
            if not d:
                continue
            tie = clave_valor(os.path.join(d, "tiempo.txt"))
            o = datos_outcar(d)
            es = energias_oszicar(d)
            runs.append({"encut": int(E), "r": int(R), "jobid": jid, "estado": tie.get("estado"),
                         "wall_s": num(tie.get("wall_s")), "n_iter": len(es), "convergido": int(o["convergido"]),
                         "e_sigma0": o["e_sigma0"], "e_oszicar": es[-1] if es else None, "npw": o["npw"],
                         "fft": "x".join(str(v) for v in o["fft"]) if o["fft"] else None,
                         "nfft": (o["fft"][0] * o["fft"][1] * o["fft"][2]) if o["fft"] else None,
                         "mem_mb": o["mem_mb"], "nions": o["nions"], "nkpts": o["nkpts"]})
    return runs


def c2_resumen(runs):
    ok = [c for c in runs if c["estado"] == "OK" and c["convergido"] and c.get("e_oszicar")]
    ref = [float(c["e_oszicar"]) for c in ok if c["encut"] == E_REF]
    if not ref:
        return [], None
    eref = ref[0]
    nions = next((c["nions"] for c in ok if c.get("nions")), None) or 23
    filas = []
    t400 = media_desv([c["wall_s"] for c in ok if c["encut"] == 400])[0]
    t250 = media_desv([c["wall_s"] for c in ok if c["encut"] == 250])[0]
    base = next((c for c in ok if c["encut"] == 250), None)
    for E in ENCUTS:
        sel = [c for c in ok if c["encut"] == E]
        if not sel:
            continue
        tm, ts = media_desv([c["wall_s"] for c in sel])
        vals = sorted(set(c["e_oszicar"] for c in sel))
        e = float(vals[0])
        err = abs(e - eref) / nions * 1000.0
        filas.append({"encut_eV": E, "n": len(sel), "t_media_s": tm, "t_desv_s": ts,
                      "t_rel_400": (tm / t400) if t400 else None, "t_rel_250": (tm / t250) if t250 else None,
                      "iteraciones": "/".join(str(c["n_iter"]) for c in sel), "npw": sel[0]["npw"],
                      "npw_rel_250": (sel[0]["npw"] / float(base["npw"])) if (base and sel[0]["npw"]) else None,
                      "fft": sel[0]["fft"], "nfft_rel_250": (sel[0]["nfft"] / float(base["nfft"])) if (base and sel[0]["nfft"]) else None,
                      "mem_mb": media_desv([c["mem_mb"] for c in sel])[0], "valores_distintos": len(vals),
                      "e_eV": vals[0], "error_meV_atomo": err if E != E_REF else None,
                      "h_nucleo": tm / 3600.0})
    return filas, eref


def c2_kpuntos():
    """convergencia con la malla de k-points: Tarea 0 (p = 1) y B.1 (mallas mas densas)"""
    filas = []
    for (K,), jid, _ in registro("jobs_t0.txt", 1):
        d = carpeta("ex_T0b_%s_%s" % (K, jid))
        if not d:
            continue
        es = energias_oszicar(d)
        tie = clave_valor(os.path.join(d, "tiempo.txt"))
        malla = re.findall(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s*$", leer(os.path.join(d, "KPOINTS")), re.M)
        o = datos_outcar(d)
        filas.append({"origen": "T0 " + K, "malla": "x".join(malla[0]) if malla else "-", "nkpts": o["nkpts"],
                      "p": 1, "e_eV": es[-1] if es else None, "wall_s": num(tie.get("wall_s"))})
    for (P, R), jid, _ in registro("jobs_B1.txt", 2):
        if R != "1":
            continue
        d = carpeta("ex_B1_p%s_r%s_%s" % (P, R, jid))
        if not d:
            continue
        es = energias_oszicar(d)
        malla = re.findall(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s*$", leer(os.path.join(d, "KPOINTS")), re.M)
        o = datos_outcar(d)
        filas.append({"origen": "B.1 p=%s" % P, "malla": "x".join(malla[0]) if malla else "-", "nkpts": o["nkpts"],
                      "p": int(P), "e_eV": es[-1] if es else None,
                      "wall_s": num(clave_valor(os.path.join(d, "tiempo.txt")).get("wall_s"))})
    ref = [f for f in filas if f["malla"] == "10x10x6" and f["e_eV"]]
    if ref:
        eref = float(ref[0]["e_eV"])
        for f in filas:
            if f["e_eV"]:
                f["error_meV_atomo"] = abs(float(f["e_eV"]) - eref) / 23.0 * 1000.0
    return filas


# ----------------------------------------------------------------------------- figura C.2
def figura_c2(plt, filas):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.edgecolor": EJE, "axes.linewidth": 0.8, "axes.labelcolor": TINTA, "xtick.color": TINTA2,
        "ytick.color": TINTA2, "text.color": TINTA, "grid.color": REJILLA, "grid.linewidth": 0.6, "axes.grid": True,
        "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "pdf.fonttype": 42, "savefig.dpi": 200})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3))
    pts = [f for f in filas if f["error_meV_atomo"] is not None and f["error_meV_atomo"] > 0]
    if not pts:                                  # sin errores positivos no hay escala logaritmica
        pts = [dict(f, error_meV_atomo=1e-3) for f in filas if f["error_meV_atomo"] is not None]
    xs = [f["t_media_s"] for f in pts]
    ys = [f["error_meV_atomo"] for f in pts]
    a1.plot(xs, ys, "-", color=C_AZUL, lw=1.2, alpha=0.6, zorder=2)
    a1.errorbar(xs, ys, xerr=[f["t_desv_s"] for f in pts], fmt="o", ms=5.5, color=C_AZUL, mec="white", mew=1.0,
                elinewidth=0.9, capsize=2, zorder=4, label=u"malla 3×3×1, 1 núcleo (3 repeticiones)")
    for f in pts:
        dx, dy, ha = 5, 4, "left"
        if f["encut_eV"] == 400:
            a1.plot([f["t_media_s"]], [f["error_meV_atomo"]], "o", ms=10, mfc="none", mec=C_NARANJA, mew=1.6, zorder=5,
                    label=u"caso base (400 eV)")
        a1.annotate("%d" % f["encut_eV"], (f["t_media_s"], f["error_meV_atomo"]), xytext=(dx, dy),
                    textcoords="offset points", fontsize=7.5, color=TINTA2, ha=ha)
    a1.axhline(1.0, ls="--", lw=0.9, color=MUTED, zorder=1)
    a1.text(xs[-1] * 1.02, 1.08, u"1 meV/átomo", fontsize=7.5, color=TINTA2, ha="right", va="bottom")
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_yticks([0.02, 0.05, 0.1, 0.2, 0.5, 1, 2])
    a1.set_yticklabels(["0.02", "0.05", "0.1", "0.2", "0.5", "1", "2"])
    a1.set_xticks([80, 100, 150, 200, 300, 500])
    a1.set_xticklabels(["80", "100", "150", "200", "300", "500"])
    a1.minorticks_off()
    a1.set_xlabel(u"tiempo de pared (s)")
    a1.set_ylabel(u"|E − E(800 eV)| por átomo (meV)")
    a1.set_title("(a)", loc="left", fontsize=9, fontweight="bold")
    a1.legend(loc="lower left")
    es = [f["encut_eV"] for f in filas]
    for clave, color, marca, et in (("t_rel_250", C_AZUL, "o", u"tiempo de pared"),
                                    ("npw_rel_250", C_NARANJA, "s", u"ondas planas por k-point"),
                                    ("nfft_rel_250", C_AQUA, "^", u"puntos de la malla FFT")):
        a2.plot(es, [f[clave] for f in filas], "-" + marca, color=color, lw=1.6, ms=5, mec="white", mew=0.8, label=et)
    a2.set_xlabel(u"ENCUT (eV)")
    a2.set_ylabel(u"crecimiento respecto de 250 eV")
    a2.set_xticks([250, 400, 550, 700, 800])
    a2.set_title("(b)", loc="left", fontsize=9, fontweight="bold")
    a2.legend(loc="upper left")
    fig.tight_layout()
    hechas = []
    os.path.isdir(FIGS) or os.makedirs(FIGS)
    for ext in ("pdf", "png"):
        ruta = os.path.join(FIGS, "C2_precision_costo." + ext)
        fig.savefig(ruta, bbox_inches="tight")
        hechas.append(ruta)
    plt.close(fig)
    return hechas


# ----------------------------------------------------------------------------- salida
def tabla_md(c1, e1, it4, red, c2, eref, kp):
    md = [u"# Parte C - precision numerica y reproducibilidad", u"",
          u"## C.1 Energia total con p procesos MPI (A.1, SCF completo, EDIFF = 1E-6 eV)", u"",
          u"E(1) = %s eV. n_eps = |E(p) - E(1)| / (eps_mach |E(1)|), con eps_mach = 2^-52 = 2.22E-16. "
          u"Resolucion de la ultima linea DAV del OSZICAR: 1E-9 eV = %.1e eps_mach |E(1)|. "
          u"EDIFF equivale a %.1e eps_mach |E(1)|." % (
              "%.9f" % e1 if e1 is not None else "-", 1e-9 / (EPS64 * abs(e1)) if e1 else float("nan"),
              EDIFF / (EPS64 * abs(e1)) if e1 else float("nan")), u"",
          u"| p | repeticiones | iteraciones SCF | valores distintos | E(p) OSZICAR (eV) | E(p) - E(1) (eV) | diferencia relativa | n_eps |",
          u"|---|---|---|---|---|---|---|---|"]
    for f in c1:
        md.append(u"| %d | %d | %s | %d | %s | %.1e | %.1e | %.1e |" % (
            f["p"], f["n"], f["iteraciones"], f["valores_distintos"], f["e_oszicar_eV"], f["de_eV"],
            f["de_relativa"], f["n_eps"]))
    md += [u"", u"Energia despues de exactamente 4 iteraciones (pasos no autoconsistentes, sin criterio de parada):", u"",
           u"| N_MPI x N_OMP | corridas | series | valores distintos | E (eV) | E - E(1x1) (eV) | n_eps |",
           u"|---|---|---|---|---|---|---|"]
    for f in it4:
        md.append(u"| %d x %d | %d | %s | %d | %s | %s | %s |" % (
            f["n_mpi"], f["n_omp"], f["corridas"], f["series"], f["valores_distintos"], f["e_iter4_eV"],
            "-" if f["de_eV"] is None else "%.1e" % f["de_eV"], "-" if f["n_eps"] is None else "%.1e" % f["n_eps"]))
    if red:
        md += [u"", u"Demostracion: suma de 2^20 numeros uniformes en [0,1) repartidos en p bloques "
                    u"(suma secuencial en cada bloque y luego de las parciales):", u"",
               u"| precision | p | diferencia con p = 1 (en eps_mach) | error frente a la suma exacta (en eps_mach) |",
               u"|---|---|---|---|"]
        for f in red:
            md.append(u"| %s | %d | %.1f | %.1f |" % (f["precision"], f["p"], f["dif_vs_p1_eps"], f["error_vs_exacta_eps"]))
    md += [u"", u"## C.2 Precision frente a costo: ENCUT (malla 3x3x1, 1 nucleo, SCF completo)", u"",
           u"Referencia: E(800 eV) = %s eV. Error por atomo = |E(ENCUT) - E(800)| / 23." % (
               "%.9f" % eref if eref is not None else "-"), u"",
           u"| ENCUT (eV) | n | t (s) | t / t(400) | iter. | ondas planas | malla FFT | memoria (MB) | E (eV) | error (meV/atomo) |",
           u"|---|---|---|---|---|---|---|---|---|---|"]
    for f in c2:
        md.append(u"| %d | %d | %.1f +/- %.1f | %.2f | %s | %s | %s | %.0f | %s | %s |" % (
            f["encut_eV"], f["n"], f["t_media_s"], f["t_desv_s"], f["t_rel_400"] or float("nan"), f["iteraciones"],
            f["npw"], f["fft"], f["mem_mb"] or float("nan"), f["e_eV"],
            "referencia" if f["error_meV_atomo"] is None else "%.3f" % f["error_meV_atomo"]))
    if kp:
        md += [u"", u"Convergencia con la malla de k-points (ENCUT = 400 eV; referencia 10x10x6 de la B.1):", u"",
               u"| origen | malla | NKPTS | p | E (eV) | error (meV/atomo) | t (s) |", u"|---|---|---|---|---|---|---|"]
        for f in kp:
            md.append(u"| %s | %s | %s | %d | %s | %s | %s |" % (
                f["origen"], f["malla"], f["nkpts"], f["p"], f["e_eV"],
                "-" if f.get("error_meV_atomo") is None else "%.4f" % f["error_meV_atomo"],
                "-" if f["wall_s"] is None else "%.0f" % f["wall_s"]))
    md.append(u"")
    with io.open(os.path.join(DATOS, "parteC_tabla.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md))


def main():
    ap = argparse.ArgumentParser(description="Analisis de la Parte C")
    ap.add_argument("--guardar", action="store_true", help="hace commit de las tablas y la figura")
    a = ap.parse_args()
    runs = c1_corridas()
    c1, e1 = c1_resumen(runs)
    it4 = c1_iteracion4(runs)
    red = demo_reduccion()
    r2 = c2_corridas()
    c2, eref = c2_resumen(r2)
    kp = c2_kpuntos()
    if not c1 or not c2:
        print("[FALLA] no encontre las corridas de la A.1 o de la C.2 (logs/ y corridas/)")
        sys.exit(1)

    print("Tarea C.1: energia total del caso base con p procesos MPI (A.1, SCF completo, EDIFF = 1E-6 eV)")
    print("E(1) = %.9f eV;  eps_mach (doble) = 2^-52 = %.3e;  eps_mach |E(1)| = %.2e eV" % (e1, EPS64, EPS64 * abs(e1)))
    print("Resolucion del OSZICAR: 1E-9 eV = %.1e eps_mach|E|;  EDIFF = 1E-6 eV = %.1e eps_mach|E|" % (
        1e-9 / (EPS64 * abs(e1)), EDIFF / (EPS64 * abs(e1))))
    print()
    print("%-4s %-6s %-9s %-8s %-22s %12s %12s %12s" % ("p", "n", "iter", "distintos", "E(p) OSZICAR (eV)", "E(p)-E(1)",
                                                         "relativa", "n_eps"))
    for f in c1:
        print("%-4d %-6d %-9s %-8d %-22s %12.1e %12.1e %12.1e" % (f["p"], f["n"], f["iteraciones"], f["valores_distintos"],
                                                                  f["e_oszicar_eV"][:22], f["de_eV"], f["de_relativa"], f["n_eps"]))
    print()
    print("Energia despues de exactamente 4 iteraciones, por N procesos MPI x T hilos OpenMP:")
    print("%-8s %-9s %-8s %-9s %-24s %12s %12s" % ("N x T", "corridas", "series", "distintos", "E (eV)", "E-E(1x1)", "n_eps"))
    for f in it4:
        print("%-8s %-9d %-8s %-9d %-24s %12s %12s" % ("%dx%d" % (f["n_mpi"], f["n_omp"]), f["corridas"], f["series"],
                                                       f["valores_distintos"], f["e_iter4_eV"][:24],
                                                       "-" if f["de_eV"] is None else "%.1e" % f["de_eV"],
                                                       "-" if f["n_eps"] is None else "%.1e" % f["n_eps"]))
    if red:
        print()
        print("Demostracion del orden de las reducciones (2^20 sumandos en p bloques):")
        for f in red:
            print("  %-17s p = %-3d  diferencia con p = 1: %8.1f eps   error frente a la suma exacta: %8.1f eps" % (
                f["precision"], f["p"], f["dif_vs_p1_eps"], f["error_vs_exacta_eps"]))
    else:
        print("(sin numpy: se omite la demostracion de la reduccion)")
    print()
    print("Tarea C.2: ENCUT (malla 3x3x1, 1 nucleo, SCF completo); referencia E(800) = %.9f eV" % eref)
    print("%-6s %-3s %16s %8s %8s %8s %-10s %8s %14s" % ("ENCUT", "n", "t (s)", "t/t400", "iter", "PW", "FFT", "mem MB",
                                                        "error meV/at"))
    for f in c2:
        print("%-6d %-3d %16s %8.2f %8s %8s %-10s %8.0f %14s" % (
            f["encut_eV"], f["n"], "%.1f +/- %.1f" % (f["t_media_s"], f["t_desv_s"]), f["t_rel_400"], f["iteraciones"],
            f["npw"], f["fft"], f["mem_mb"] or 0, "referencia" if f["error_meV_atomo"] is None else "%.3f" % f["error_meV_atomo"]))
    distintos = [f for f in c2 if f["valores_distintos"] != 1]
    print("Las 3 repeticiones de cada ENCUT dan %s." % (
        "la misma energia (12 cifras)" if not distintos else "energias distintas en %d casos" % len(distintos)))
    if kp:
        print()
        print("Malla de k-points (ENCUT = 400 eV; referencia 10x10x6 de la B.1):")
        for f in kp:
            print("  %-9s %-8s NKPTS = %-3s p = %-3d E = %-22s error = %s meV/atomo" % (
                f["origen"], f["malla"], f["nkpts"], f["p"], f["e_eV"],
                "-" if f.get("error_meV_atomo") is None else "%.4f" % f["error_meV_atomo"]))

    escribir_csv("parteC1_energias.csv", ["serie", "n_mpi", "n_omp", "r", "jobid", "estado", "n_iter", "convergido",
                                          "e_final_oszicar", "e_iter4_oszicar", "e_sigma0_outcar", "ranks", "hilos"], runs)
    escribir_csv("parteC1_resumen.csv", ["p", "n", "iteraciones", "valores_distintos", "e_oszicar_eV", "de_eV",
                                         "de_relativa", "n_eps", "e_outcar_eV", "de_outcar_eV", "n_eps_outcar"], c1)
    escribir_csv("parteC1_iteracion4.csv", ["n_mpi", "n_omp", "corridas", "series", "valores_distintos", "e_iter4_eV",
                                            "de_eV", "n_eps", "bandas_por_proceso"], it4)
    if red:
        escribir_csv("parteC1_reduccion.csv", ["precision", "eps_mach", "n_terminos", "p", "suma", "dif_vs_p1_rel",
                                               "dif_vs_p1_eps", "error_vs_exacta_eps"], red)
    escribir_csv("parteC2_corridas.csv", ["encut", "r", "jobid", "estado", "wall_s", "n_iter", "convergido", "e_sigma0",
                                          "e_oszicar", "npw", "fft", "nfft", "mem_mb", "nions", "nkpts"], r2)
    escribir_csv("parteC2_resumen.csv", ["encut_eV", "n", "t_media_s", "t_desv_s", "t_rel_400", "t_rel_250", "iteraciones",
                                         "npw", "npw_rel_250", "fft", "nfft_rel_250", "mem_mb", "valores_distintos", "e_eV",
                                         "error_meV_atomo", "h_nucleo"], c2)
    escribir_csv("parteC2_kpuntos.csv", ["origen", "malla", "nkpts", "p", "e_eV", "error_meV_atomo", "wall_s"], kp)
    tabla_md(c1, e1, it4, red, c2, eref, kp)
    print()
    print("Tablas guardadas en datos/parteC1_*.csv, datos/parteC2_*.csv y datos/parteC_tabla.md")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[AVISO] falta matplotlib: no se hizo la figura (python3 -m pip install --user matplotlib)")
        sys.exit(1)
    hechas = figura_c2(plt, c2)
    print("Figura guardada: " + ", ".join(os.path.relpath(h, REPO) for h in hechas))
    if a.guardar:
        out = git("add", "datos", "reporte/figuras", "analisis/parteC_analisis.py")
        if out.startswith("ERROR"):
            print("AVISO: git add fallo: " + out)
            return
        if git("diff", "--cached", "--quiet").startswith("ERROR"):
            out = git("commit", "-q", "-m", "Parte C: reproducibilidad con p (C.1) y precision frente a costo con ENCUT (C.2)")
            print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
        else:
            print("Nada nuevo que guardar en Git.")


if __name__ == "__main__":
    main()
