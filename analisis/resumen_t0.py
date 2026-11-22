#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
resumen_t0.py - resume la calibracion de la Tarea 0 y fija el caso base.

Uso (desde la carpeta del repositorio):
  python3 analisis/resumen_t0.py                 # tabla, verificaciones y recomendacion
  python3 analisis/resumen_t0.py --fijar K5      # fija K5 como caso base + tabla de la Tarea 0 + commit
  python3 analisis/resumen_t0.py --nkpts 8 8 2   # k-points irreducibles esperados de una malla

Regla de eleccion del caso base: entre los candidatos que terminaron bien y cuyo
t1 (1 nucleo) esta entre 15 y 60 min, el mas cercano a 30 min (en escala
logaritmica, el centro de la ventana). Asi queda margen de un factor 2 hacia
ambos lados por la variabilidad del cluster.

Solo usa la biblioteca estandar de Python. Todo lo que imprime es ASCII.
"""
from __future__ import print_function
import argparse
import csv
import io
import math
import os
import re
import shutil
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REG = os.path.join(REPO, "logs", "jobs_t0.txt")
NELECT_ESPERADO = 202.0      # 8 W_sv (14 e) + 15 Se (6 e)
NIONS_ESPERADO = 23
T_MIN, T_MAX, T_OBJ = 15.0, 60.0, 30.0   # minutos

# k-points irreducibles de mallas Gamma n x n x m para ESTE POSCAR (simetria C3v de la
# vacancia + inversion temporal), calculados con el POSCAR generado. VASP debe coincidir.
TABLA_NKPTS = {
    (2, 2, 1): 2, (2, 2, 2): 4, (2, 2, 3): 4, (2, 2, 4): 6, (3, 3, 1): 3, (3, 3, 2): 6,
    (3, 3, 3): 7, (3, 3, 4): 10, (4, 4, 1): 4, (4, 4, 2): 8, (4, 4, 3): 9, (4, 4, 4): 13,
    (5, 5, 1): 5, (5, 5, 2): 10, (5, 5, 3): 12, (5, 5, 4): 17, (6, 6, 1): 7, (6, 6, 2): 14,
    (6, 6, 3): 17, (6, 6, 4): 24, (7, 7, 1): 8, (7, 7, 2): 16, (7, 7, 3): 20, (7, 7, 4): 28,
    (8, 8, 1): 10, (8, 8, 2): 20, (8, 8, 3): 25, (8, 8, 4): 35, (9, 9, 1): 12, (9, 9, 2): 24,
    (9, 9, 3): 31, (9, 9, 4): 43, (10, 10, 1): 14, (10, 10, 2): 28, (10, 10, 3): 36,
    (10, 10, 4): 50, (11, 11, 1): 16, (11, 11, 2): 32, (11, 11, 3): 42, (11, 11, 4): 58,
    (12, 12, 1): 19, (12, 12, 2): 38, (12, 12, 3): 50, (12, 12, 4): 69, (13, 13, 1): 21,
    (13, 13, 2): 42, (13, 13, 3): 56, (13, 13, 4): 77, (14, 14, 1): 24, (14, 14, 2): 48,
    (14, 14, 3): 64, (14, 14, 4): 88,
}


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


def num(x, tipo=float):
    try:
        return tipo(x)
    except (TypeError, ValueError):
        return None


def buscar(patron, texto, grupo=1, ultimo=False, tipo=str):
    m = re.findall(patron, texto)
    if not m:
        return None
    v = m[-1] if ultimo else m[0]
    if isinstance(v, tuple):
        v = v[grupo - 1]
    return num(v, tipo) if tipo is not str else v


def parametro_incar(texto, tag):
    for l in texto.splitlines():
        l = l.split("#")[0].split("!")[0]
        m = re.match(r"\s*" + tag + r"\s*=\s*(\S+)", l, re.I)
        if m:
            return m.group(1)
    return None


def git(*args):
    try:
        return subprocess.check_output(("git",) + args, cwd=REPO, stderr=subprocess.STDOUT).decode().strip()
    except (OSError, subprocess.CalledProcessError) as e:
        out = getattr(e, "output", b"")
        return "ERROR: " + (out.decode().strip() if out else str(e))


def estado_slurm(jobid):
    """estado del job segun Slurm: PENDING/RUNNING si esta en cola; si no, el de sacct (TIMEOUT, FAILED...)"""
    for cmd in (["squeue", "-h", "-j", str(jobid), "-o", "%T"], ["sacct", "-n", "-X", "-j", str(jobid), "-o", "State%20"]):
        try:
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT).decode().split()
        except (OSError, subprocess.CalledProcessError):
            continue
        if out:
            return out[0]
    return None


def sacct_maxrss_mb(jobid):
    """MaxRSS mas grande de los pasos del job, en MB (None si sacct no responde)."""
    try:
        out = subprocess.check_output(["sacct", "-j", str(jobid), "-n", "-P", "-o", "JobID,MaxRSS"],
                                      stderr=subprocess.STDOUT).decode()
    except (OSError, subprocess.CalledProcessError):
        return None
    mejor = None
    fac = {"K": 1.0 / 1024, "M": 1.0, "G": 1024.0, "T": 1024.0 ** 2}
    for l in out.splitlines():
        p = l.split("|")
        if len(p) < 2 or not p[1].strip():
            continue
        m = re.match(r"([\d.]+)([KMGT]?)", p[1].strip())
        if m:
            v = float(m.group(1)) * fac.get(m.group(2), 1.0 / 1024 ** 2)
            mejor = v if mejor is None else max(mejor, v)
    return mejor


# ----------------------------------------------------------------------------- lectura
def candidatos():
    c = {}
    orden = []
    for l in leer(os.path.join(REPO, "entradas", "candidatos.txt")).splitlines():
        p = l.split()
        if not p or p[0].startswith("#") or len(p) < 4:
            continue
        malla = (int(p[1]), int(p[2]), int(p[3]))
        nk = num(p[4], int) if len(p) > 4 else None
        if nk is None:
            nk = TABLA_NKPTS.get(malla)
        c[p[0]] = {"malla": malla, "nk_esperado": nk, "equiv": p[5] if len(p) > 5 else "?"}
        orden.append(p[0])
    return c, orden


def leer_outcar(ruta):
    t = leer(ruta)
    o = {"existe": bool(t)}
    if not t:
        return o
    o["version"] = buscar(r"(vasp\.[0-9.]+[^\n]*)", t)
    o["ranks"] = buscar(r"running\s+(\d+)\s+mpi-ranks", t, tipo=int)
    o["hilos"] = buscar(r"with\s+(\d+)\s+threads/rank", t, tipo=int)
    if o["ranks"] is None:
        o["ranks"] = buscar(r"running on\s+(\d+)\s+total cores", t, tipo=int)
    o["nions"] = buscar(r"NIONS\s*=\s*(\d+)", t, tipo=int)
    o["nkpts"] = buscar(r"NKPTS\s*=\s*(\d+)", t, tipo=int)
    # la linea "number of bands NBANDS=" es el valor que VASP usa de verdad (puede redondearlo)
    o["nbands"] = buscar(r"number of bands\s+NBANDS\s*=\s*(\d+)", t, ultimo=True, tipo=int)
    if o["nbands"] is None:
        o["nbands"] = buscar(r"NBANDS\s*=\s*(\d+)", t, ultimo=True, tipo=int)
    o["nelect"] = buscar(r"NELECT\s*=\s*([\d.]+)", t, tipo=float)
    o["encut"] = buscar(r"ENCUT\s*=\s*([\d.]+)\s*eV", t, tipo=float)
    g = re.search(r"dimension x,y,z NGX =\s*(\d+)\s+NGY =\s*(\d+)\s+NGZ =\s*(\d+)", t)
    o["fft"] = "%s x %s x %s" % g.groups() if g else None
    pw = re.search(r"maximum and minimum number of plane-waves per node\s*:\s*(\d+)\s+(\d+)", t)
    o["pw_max"] = int(pw.group(1)) if pw else buscar(r"maximum number of plane-waves\s*:\s*(\d+)", t, tipo=int)
    o["mem_est_mb"] = None
    v = buscar(r"total amount of memory used by VASP MPI-rank0\s+([\d.]+)\s*kBytes", t, tipo=float)
    if v is not None:
        o["mem_est_mb"] = v / 1024.0
    v = buscar(r"Maximum memory used \(kb\):\s*([\d.]+)", t, tipo=float)
    o["mem_max_mb"] = v / 1024.0 if v is not None else None
    o["elapsed_s"] = buscar(r"Elapsed time \(sec\):\s*([\d.]+)", t, tipo=float)
    loops = [float(x) for x in re.findall(r"LOOP:\s+cpu time\s+[\d.]+\s*:\s*real time\s+([\d.]+)", t)]
    o["n_iter"] = len(loops)
    o["t_iter_s"] = sum(loops) / len(loops) if loops else None
    o["convergido"] = "aborting loop because EDIFF is reached" in t
    o["fin_normal"] = "General timing and accounting" in t
    o["e_sigma0"] = buscar(r"energy\(sigma->0\)\s*=\s*([-\d.]+)", t, ultimo=True, tipo=float)
    o["toten"] = buscar(r"free\s+energy\s+TOTEN\s*=\s*([-\d.]+)\s*eV", t, ultimo=True, tipo=float)
    return o


def leer_carga(ruta):
    """maximo de CPUs ocupados en el socket del job durante la corrida ('x/y')"""
    mejor, total = None, None
    for l in leer(ruta).splitlines():
        m = re.search(r"cpus_ocupadas_socket=(\d+)/(\d+)", l)
        if m:
            v = int(m.group(1))
            total = int(m.group(2))
            mejor = v if mejor is None else max(mejor, v)
    return None if mejor is None else "%d/%d" % (mejor, total)


def corridas():
    """lee logs/jobs_t0.txt y devuelve, por candidato, la ultima corrida (prefiere las OK)"""
    if not os.path.isfile(REG):
        return {}
    porcand = {}
    for l in leer(REG).splitlines():
        p = l.split()
        if len(p) >= 2:
            porcand.setdefault(p[0], []).append(p[1])
    res = {}
    for k, ids in porcand.items():
        elegido = None
        for jid in reversed(ids):
            d = os.path.join(REPO, "corridas", "ex_T0b_%s_%s" % (k, jid))
            tie = clave_valor(os.path.join(d, "tiempo.txt"))
            r = {"cand": k, "jobid": jid, "dir": d, "existe": os.path.isdir(d), "tiempo": tie,
                 "entorno": clave_valor(os.path.join(d, "entorno.txt"))}
            if elegido is None:
                elegido = r
            if tie.get("estado") == "OK":
                elegido = r
                break
        res[k] = elegido
    return res


# ----------------------------------------------------------------------------- analisis
def analizar(incar):
    cand, orden = candidatos()
    runs = corridas()
    encut_in = num(parametro_incar(incar, "ENCUT"), float)
    nbands_in = num(parametro_incar(incar, "NBANDS"), int)
    filas, avisos = [], []
    for k in orden:
        c = cand[k]
        r = runs.get(k)
        f = {"cand": k, "malla": "%dx%dx%d" % c["malla"], "equiv": c["equiv"], "nk_esp": c["nk_esperado"],
             "jobid": None, "estado": "SIN_ENVIAR", "ok": False}
        if r is None:
            filas.append(f)
            continue
        f["jobid"] = r["jobid"]
        if not r["existe"]:
            st = estado_slurm(r["jobid"])
            f["estado"] = "EN_COLA" if st in (None, "PENDING", "CONFIGURING") else st
            filas.append(f)
            continue
        o = leer_outcar(os.path.join(r["dir"], "OUTCAR"))
        tie = r["tiempo"]
        f.update(o)
        f["nodo"] = r["entorno"].get("nodo")
        f["git_hash"] = r["entorno"].get("git_hash")
        f["carga"] = leer_carga(os.path.join(r["dir"], "carga.txt"))
        if tie.get("estado"):
            f["estado"] = tie["estado"]
        else:
            st = estado_slurm(r["jobid"])
            f["estado"] = "CORRIENDO" if st in (None, "RUNNING", "COMPLETING") else st
        wall = num(tie.get("wall_s"))
        f["wall_s"] = wall
        cand_t = [x for x in (wall, o.get("elapsed_s")) if x]
        f["t1_min"] = max(cand_t) / 60.0 if cand_t else None
        if o.get("t_iter_s") and o.get("nkpts"):
            f["t_k_s"] = o["t_iter_s"] / o["nkpts"]
        f["ok"] = f["estado"] == "OK" and bool(o.get("convergido")) and f["t1_min"] is not None
        # ---- verificaciones contra lo esperado
        if f["estado"] == "OK":
            chk = [("NIONS", o.get("nions"), NIONS_ESPERADO), ("NELECT", o.get("nelect"), NELECT_ESPERADO),
                   ("NBANDS", o.get("nbands"), nbands_in), ("ENCUT", o.get("encut"), encut_in),
                   ("procesos MPI", o.get("ranks"), 1), ("hilos OpenMP", o.get("hilos") or 1, 1)]
            for nombre, v, esp in chk:
                if v is None or esp is None:
                    avisos.append("%s: no pude leer %s del OUTCAR" % (k, nombre))
                elif abs(float(v) - float(esp)) > 1e-6:
                    avisos.append("%s: %s = %s y se esperaba %s (no se puede usar como caso base)" % (k, nombre, v, esp))
                    f["ok"] = False
            nk = o.get("nkpts")
            if c["nk_esperado"] and nk is not None and nk != c["nk_esperado"]:
                avisos.append("%s: VASP uso %s k-points irreducibles y se esperaban %s (simetria distinta a la "
                              "prevista; mandame una captura antes de fijar)" % (k, nk, c["nk_esperado"]))
        filas.append(f)
    return filas, avisos


def elegir(filas):
    dentro = [f for f in filas if f["ok"] and T_MIN <= f["t1_min"] <= T_MAX]
    if not dentro:
        return None
    return sorted(dentro, key=lambda f: (abs(math.log(f["t1_min"] / T_OBJ)), f.get("nkpts") or 0))[0]


def predecir(filas):
    """Si ningun candidato quedo en la ventana: ajusta t_iter = a + b*NKPTS y recomienda una malla."""
    buenos = [f for f in filas if f["ok"] and f.get("t_iter_s") and f.get("nkpts")]
    if not buenos:
        return None
    xs = [f["nkpts"] for f in buenos]
    ys = [f["t_iter_s"] for f in buenos]
    if len(set(xs)) >= 2:
        n = float(len(xs))
        mx, my = sum(xs) / n, sum(ys) / n
        sxx = sum((x - mx) ** 2 for x in xs)
        b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
        a = my - b * mx
        if b <= 0:
            a, b = 0.0, my / mx
    else:
        a, b = 0.0, ys[0] / xs[0]
    its = sorted(f["n_iter"] for f in buenos)
    nit = its[len(its) // 2]
    mejor = None
    for malla, nk in sorted(TABLA_NKPTS.items()):
        t = nit * (a + b * nk) / 60.0
        d = abs(math.log(t / T_OBJ)) if t > 0 else 1e9
        if mejor is None or d < mejor[0]:
            mejor = (d, malla, nk, t)
    return {"a": a, "b": b, "nit": nit, "malla": mejor[1], "nk": mejor[2], "t1": mejor[3]}


def fmt(v, f="%.1f", nada="-"):
    return nada if v is None else (f % v)


def imprimir(filas, avisos):
    print("Calibracion de la Tarea 0: WSe2 bulk 2x2x1 con V_Se, SCF en 1 nucleo")
    print()
    cab = "%-4s %-8s %-8s %5s %4s %8s %8s %8s %8s %7s  %s" % (
        "Cand", "Malla", "Job", "NKPTS", "Iter", "t/iter", "t/k-pt", "t1[min]", "Mem[MB]", "Socket", "Estado")
    print(cab)
    print("-" * len(cab))
    for f in filas:
        nota = f["estado"]
        if f["estado"] == "OK" and not f["ok"]:
            nota = "REVISAR (ver AVISO abajo)"
        if f["ok"]:
            if f["t1_min"] < T_MIN:
                nota = "OK (< 15 min)"
            elif f["t1_min"] > T_MAX:
                nota = "OK (> 60 min)"
            else:
                nota = "OK (en la ventana)"
        print("%-4s %-8s %-8s %5s %4s %8s %8s %8s %8s %7s  %s" % (
            f["cand"], f["malla"], f["jobid"] or "-", fmt(f.get("nkpts"), "%d"), fmt(f.get("n_iter"), "%d"),
            fmt(f.get("t_iter_s"), "%.1f s"), fmt(f.get("t_k_s"), "%.2f s"), fmt(f.get("t1_min"), "%.1f"),
            fmt(f.get("mem_max_mb"), "%.0f"), f.get("carga") or "-", nota))
    print()
    print("Socket = maximo de CPUs logicos ocupados (>50 %) en el socket del job durante la corrida,")
    print("         incluido el propio job. Valores altos indican competencia con otros jobs.")
    if avisos:
        print()
        for a in avisos:
            print("AVISO: " + a)
    print()


def recomendar(filas):
    e = elegir(filas)
    terminados = [f for f in filas if f["estado"] not in ("SIN_ENVIAR", "EN_COLA", "CORRIENDO")]
    pendientes = [f for f in filas if f["estado"] in ("EN_COLA", "CORRIENDO")]
    if pendientes:
        print("Todavia hay candidatos en cola o corriendo: %s" % ", ".join(f["cand"] for f in pendientes))
    if e:
        print("RECOMENDACION: %s (malla %s, %d k-points): t1 = %.1f min, el mas cercano a 30 min"
              % (e["cand"], e["malla"], e["nkpts"], e["t1_min"]))
        print("               dentro de la ventana de 15 a 60 min.")
        print("Para fijarlo:  python3 analisis/resumen_t0.py --fijar %s" % e["cand"])
    elif terminados:
        p = predecir(filas)
        if p:
            k = "K%d" % (max([int(f["cand"][1:]) for f in filas if f["cand"][1:].isdigit()] + [0]) + 1)
            print("Ningun candidato quedo entre 15 y 60 min.")
            print("Con los tiempos medidos (%.2f s por k-point e iteracion, %d iteraciones) la malla"
                  % (p["b"], p["nit"]))
            print("%dx%dx%d (%d k-points) deberia tardar unos %.0f min. Para probarla:"
                  % (p["malla"] + (p["nk"], p["t1"])))
            print("  bash slurm/nuevo_candidato.sh %s %d %d %d" % ((k,) + p["malla"]))
        else:
            print("Ningun candidato termino bien: revisa los logs/ y mandame una captura.")
    return e


def guardar_csv(filas):
    ruta = os.path.join(REPO, "datos", "tarea0_calibracion.csv")
    if not os.path.isdir(os.path.dirname(ruta)):
        os.makedirs(os.path.dirname(ruta))
    cols = ["cand", "malla", "equiv", "jobid", "nodo", "estado", "nkpts", "nk_esp", "n_iter", "t_iter_s", "t_k_s",
            "wall_s", "elapsed_s", "t1_min", "mem_est_mb", "mem_max_mb", "pw_max", "fft", "nions", "nelect",
            "nbands", "encut", "e_sigma0", "toten", "carga", "version"]
    with io.open(ruta, "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for f in filas:
            w.writerow(["" if f.get(c) is None else f.get(c) for c in cols])
    return ruta


# ----------------------------------------------------------------------------- fijar
def fijar(filas, nombre, incar, forzar=False):
    f = [x for x in filas if x["cand"] == nombre]
    if not f:
        sys.exit("ERROR: no existe el candidato %s" % nombre)
    f = f[0]
    if not f["ok"]:
        sys.exit("ERROR: %s no termino bien (estado %s); no se puede fijar" % (nombre, f["estado"]))
    if not (T_MIN <= f["t1_min"] <= T_MAX) and not forzar:
        sys.exit("ERROR: t1 de %s = %.1f min esta fuera de 15-60 min (usa --forzar si el docente lo autoriza)"
                 % (nombre, f["t1_min"]))
    run = os.path.join(REPO, "corridas", "ex_T0b_%s_%s" % (nombre, f["jobid"]))
    src = os.path.join(REPO, "entradas", "candidatos", nombre)
    dst = os.path.join(REPO, "entradas", "caso_base")
    if not os.path.isdir(dst):
        os.makedirs(dst)
    for a in ("INCAR", "KPOINTS", "POSCAR"):
        shutil.copy(os.path.join(src, a), os.path.join(dst, a))
    with io.open(os.path.join(dst, "ORIGEN.txt"), "w", encoding="utf-8") as fh:
        fh.write(u"candidato = %s\nmalla = %s\njob = %s\nt1_min = %.2f\nnodo = %s\ncommit_entradas = %s\n"
                 % (nombre, f["malla"], f["jobid"], f["t1_min"], f.get("nodo"), f.get("git_hash")))
    ev = os.path.join(REPO, "datos", "tarea0", nombre)
    if not os.path.isdir(ev):
        os.makedirs(ev)
    for a in ("OSZICAR", "tiempo.txt", "entorno.txt", "modulos.txt", "ldd.txt", "entradas.md5",
              "carga.txt", "variables.txt"):
        if os.path.isfile(os.path.join(run, a)):
            shutil.copy(os.path.join(run, a), os.path.join(ev, a))
    lineas = leer(os.path.join(run, "OUTCAR")).splitlines()
    with io.open(os.path.join(ev, "OUTCAR_extracto.txt"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(lineas[:400]) + u"\n\n[... recorte ...]\n\n" + u"\n".join(lineas[-60:]) + u"\n")

    geo = clave_valor(os.path.join(REPO, "entradas", "geometria_bulk_2x2x1_VSe.txt"))
    pot = clave_valor(os.path.join(REPO, "entradas", "POTCAR.titel"))
    rss = sacct_maxrss_mb(f["jobid"])
    mods = re.findall(r"\d+\)\s+(\S+)", leer(os.path.join(run, "modulos.txt")))
    mod = "modulos " + ", ".join(mods) if mods else leer(os.path.join(run, "modulos.txt")).strip().replace("\n", " ")
    filas_t = [
        ("Sistema fisico", "2H-WSe2 bulk, supercelda 2x2x1 con una vacancia de Se (V_Se), el defecto de mi tesis"),
        ("Geometria", "experimental (Schutte et al. 1987): a = %s A, c = %s A, z(Se) = %s; sin relajar"
         % (geo.get("a_primitiva_A", "?"), geo.get("c_A", "?"), geo.get("z_Se", "?"))),
        ("Numero de atomos", "%s (%s W + %s Se)" % (f.get("nions"), geo.get("n_W", "?"), geo.get("n_Se", "?"))),
        ("Electrones de valencia (NELECT)", fmt(f.get("nelect"), "%.0f")),
        ("PAW / funcional", "%s; %s / PBE (GGA = PE)" % (pot.get("TITEL 1", "?"), pot.get("TITEL 2", "?"))),
        ("Calculo", "SCF estatico desde cero (ISTART = 0, ICHARG = 2, IBRION = -1, NSW = 0), EDIFF = %s, "
                    "ISMEAR = %s, SIGMA = %s, ISPIN = %s" % (parametro_incar(incar, "EDIFF"),
                    parametro_incar(incar, "ISMEAR"), parametro_incar(incar, "SIGMA"), parametro_incar(incar, "ISPIN"))),
        ("ENCUT / PREC / LREAL", "%s eV / %s / %s" % (fmt(f.get("encut"), "%.0f"), parametro_incar(incar, "PREC"),
                                                     parametro_incar(incar, "LREAL"))),
        ("Malla de k-points", "Gamma %s (= %s en la primitiva) -> %s k-points irreducibles"
         % (f["malla"], f["equiv"], f.get("nkpts"))),
        ("NBANDS", "%s (fijo, multiplo de 48: no cambia con el numero de procesos)" % f.get("nbands")),
        ("Ondas planas por k-point (max)", fmt(f.get("pw_max"), "%d")),
        ("Malla FFT (NGX x NGY x NGZ)", f.get("fft") or "-"),
        ("Iteraciones SCF", fmt(f.get("n_iter"), "%d")),
        ("t(1 nucleo)", "%.1f min (pared %s s; Elapsed de VASP %s s)" % (f["t1_min"], fmt(f.get("wall_s"), "%.0f"),
                                                                     fmt(f.get("elapsed_s"), "%.0f"))),
        ("Memoria estimada por VASP", "%s MB" % fmt(f.get("mem_est_mb"), "%.0f")),
        ("Memoria medida", "%s MB (OUTCAR) / %s MB (sacct MaxRSS)" % (fmt(f.get("mem_max_mb"), "%.0f"),
                                                                     fmt(rss, "%.0f"))),
        ("Limite de memoria de la QOS", "62 600 MB"),
        ("Software", "%s; %s; binario vasp_std" % ((f.get("version") or "?").strip(), mod or "?")),
        ("Nodo / CPUs ocupados en su socket", "%s / %s" % (f.get("nodo"), f.get("carga") or "-")),
        ("Semilla aleatoria", "No aplica: VASP inicia las funciones de onda con INIWAV = 1 y no expone una semilla "
                              "(RANDOM_SEED solo afecta la dinamica molecular). La reproducibilidad se mide en la C.1"),
        ("Job / commit de las entradas", "%s / %s" % (f["jobid"], f.get("git_hash"))),
    ]
    md = [u"# Tarea 0 - caso base", u"", u"| Propiedad | Valor |", u"|---|---|"]
    md += [u"| %s | %s |" % (a, b) for a, b in filas_t]
    md += [u"", u"Calibracion completa: `datos/tarea0_calibracion.csv`. Evidencia de la corrida: `datos/tarea0/%s/`."
           % nombre, u""]
    ruta_md = os.path.join(REPO, "datos", "tarea0_tabla_caso_base.md")
    with io.open(ruta_md, "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md))
    print("Caso base fijado: entradas/caso_base/ (INCAR, KPOINTS, POSCAR, ORIGEN.txt)")
    print("Tabla de la Tarea 0: datos/tarea0_tabla_caso_base.md")
    print()
    for a, b in filas_t:
        print("  %-34s %s" % (a, b))
    print()
    out = git("add", "entradas/caso_base", "datos")
    if out.startswith("ERROR"):
        print("AVISO: git add fallo: " + out)
        return
    if git("diff", "--cached", "--quiet").startswith("ERROR"):
        out = git("commit", "-q", "-m", "Tarea 0: caso base fijado (%s, malla %s, t1 = %.1f min)"
                  % (nombre, f["malla"], f["t1_min"]))
        if out.startswith("ERROR"):
            print("AVISO: git commit fallo: " + out)
        else:
            print("Commit hecho: " + git("log", "--oneline", "-1"))
    else:
        print("Nada nuevo que guardar en Git.")


# ----------------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description="Resumen de la calibracion de la Tarea 0")
    p.add_argument("--fijar", metavar="KX", help="fija ese candidato como caso base")
    p.add_argument("--forzar", action="store_true", help="permite fijar un candidato fuera de 15-60 min")
    p.add_argument("--nkpts", nargs=3, type=int, metavar=("N1", "N2", "N3"),
                   help="imprime los k-points irreducibles esperados de una malla y termina")
    a = p.parse_args()
    if a.nkpts:
        v = TABLA_NKPTS.get(tuple(a.nkpts))
        print(v if v is not None else "?")
        return
    incar = leer(os.path.join(REPO, "entradas", "INCAR_base"))
    filas, avisos = analizar(incar)
    if not any(f["jobid"] for f in filas):
        print("Todavia no hay corridas: primero  bash slurm/enviar.sh")
        return
    imprimir(filas, avisos)
    ruta = guardar_csv(filas)
    print("Tabla guardada en " + os.path.relpath(ruta, REPO))
    print()
    if a.fijar:
        fijar(filas, a.fijar, incar, a.forzar)
    else:
        recomendar(filas)


if __name__ == "__main__":
    main()
