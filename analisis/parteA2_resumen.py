#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteA2_resumen.py - Tarea A.2: escalamiento fuerte con OpenMP (memoria compartida).

Lee las corridas de logs/jobs_A2.txt: 1 proceso MPI con p hilos, p = 1..24, 3
repeticiones, con una carga fija de 4 iteraciones SCF (NELM = 4) del caso base.
Verifica cada una y reporta:
  * t(p) de la corrida (media y desviacion estandar), S(p) = t(1)/t(p), E(p) = S/p,
    con t(1) medido con la misma carga (p = 1 de este barrido);
  * t por iteracion SCF (LOOP del OUTCAR) y su S y E;
  * la comparacion por iteracion con MPI: las primeras 4 iteraciones de las corridas
    de la A.1 (mismo caso, mismo trabajo por iteracion).

Uso (desde la carpeta del repositorio):
  python3 analisis/parteA2_resumen.py             # tablas, verificaciones y CSV
  python3 analisis/parteA2_resumen.py --guardar   # ademas copia la evidencia a datos/ y hace commit

Archivos que escribe:
  datos/parteA2_corridas.csv     una fila por corrida
  datos/parteA2_resumen.csv      una fila por p (t, S, E de la corrida y por iteracion, memoria)
  datos/parteA2_mpi_vs_omp.csv   S y E por iteracion: MPI (A.1) y OpenMP (A.2)
  datos/parteA2_tabla.md         tablas para el reporte

Una corrida se puede excluir agregando su numero de job y el motivo a
datos/parteA2_excluidas.txt; luego  bash slurm/parteA2_enviar.sh --faltantes  la repite.

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

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P_LISTA = [1, 2, 4, 8, 12, 16, 24]
NELM_FIJO = 4
CARGA_ALTA = 0.25
DESVIO = 0.15
DESCRIPCION = ("1 proceso con p hilos OpenMP (OMP_PLACES=cores, OMP_PROC_BIND=close), un nucleo fisico por hilo; "
               "nucleos pedidos como en la A.1 y paso de srun elegido por su mascara real; OPENBLAS_NUM_THREADS=1; "
               "carga fija de 4 iteraciones SCF (NELM = 4)")


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


def buscar(patron, texto, ultimo=False, tipo=str):
    m = re.findall(patron, texto)
    if not m:
        return None
    v = m[-1] if ultimo else m[0]
    if isinstance(v, tuple):
        v = v[0]
    return num(v, tipo) if tipo is not str else v


def git(*args):
    try:
        return subprocess.check_output(("git",) + args, cwd=REPO, stderr=subprocess.STDOUT).decode().strip()
    except (OSError, subprocess.CalledProcessError) as e:
        out = getattr(e, "output", b"")
        return "ERROR: " + (out.decode().strip() if out else str(e))


def estado_slurm(jobid):
    for cmd in (["squeue", "-h", "-j", str(jobid), "-o", "%T"], ["sacct", "-n", "-X", "-j", str(jobid), "-o", "State%20"]):
        try:
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT).decode().split()
        except (OSError, subprocess.CalledProcessError):
            continue
        if out:
            return out[0]
    return None


def media_desv(xs):
    n = len(xs)
    if n == 0:
        return None, None
    m = sum(xs) / float(n)
    if n == 1:
        return m, 0.0
    return m, math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))


def mediana(xs):
    s = sorted(xs)
    n = len(s)
    return None if n == 0 else (s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2]))


def fmt(v, f="%.1f", nada="-"):
    return nada if v is None else (f % v)


def cociente(a, sa, b, sb):
    """a/b con su incertidumbre propagada (desviaciones relativas en cuadratura)"""
    if not a or not b:
        return None, None
    q = a / b
    return q, q * math.sqrt(((sa or 0.0) / a) ** 2 + ((sb or 0.0) / b) ** 2)


# ----------------------------------------------------------------------------- lectura
def leer_outcar(ruta):
    t = leer(ruta)
    o = {}
    if not t:
        return o
    o["version"] = buscar(r"(vasp\.[0-9.]+[^\n]*)", t)
    m = re.search(r"running\s+(\d+)\s+mpi-ranks,\s+with\s+(\d+)\s+threads/rank", t)
    o["ranks"] = int(m.group(1)) if m else None
    o["hilos_vasp"] = int(m.group(2)) if m else None
    o["nkpts"] = buscar(r"NKPTS\s*=\s*(\d+)", t, tipo=int)
    o["nbands"] = buscar(r"number of bands\s+NBANDS\s*=\s*(\d+)", t, ultimo=True, tipo=int)
    pares = [(float(a), float(b)) for a, b in re.findall(r"LOOP:\s+cpu time\s+([\d.]+)\s*:\s*real time\s+([\d.]+)", t)]
    o["loops"] = [b for a, b in pares]
    o["n_iter"] = len(pares)
    o["t_iter_s"] = sum(b for a, b in pares) / len(pares) if pares else None
    o["cpu_real"] = (sum(a for a, b in pares) / sum(b for a, b in pares)) if pares and sum(b for a, b in pares) > 0 else None
    o["elapsed_s"] = buscar(r"Elapsed time \(sec\):\s*([\d.]+)", t, tipo=float)
    v = buscar(r"Maximum memory used \(kb\):\s*([\d.]+)", t, tipo=float)
    o["mem_mb"] = v / 1024.0 if v is not None else None
    o["convergido"] = "aborting loop because EDIFF is reached" in t
    o["toten"] = buscar(r"free\s+energy\s+TOTEN\s*=\s*([-\d.]+)\s*eV", t, ultimo=True)
    o["e_sigma0"] = buscar(r"energy\(sigma->0\)\s*=\s*([-\d.]+)", t, ultimo=True)
    return o


def excluidas(ruta):
    d = {}
    for l in leer(ruta).splitlines():
        p = l.split(None, 1)
        if p and not p[0].startswith("#"):
            d[p[0]] = p[1].strip() if len(p) > 1 else "sin motivo"
    return d


def corridas(reg, prefijo, excl, origen, nelm_fijo=None):
    if not os.path.isfile(reg):
        return []
    exc = excluidas(excl)
    res = []
    for orden, l in enumerate(leer(reg).splitlines(), 1):
        p = l.split()
        if len(p) < 3 or not p[0].isdigit():
            continue
        P, R, jid = int(p[0]), int(p[1]), p[2]
        d = os.path.join(REPO, "corridas", "%s_p%d_r%d_%s" % (prefijo, P, R, jid))
        c = {"orden": orden, "p": P, "r": R, "jobid": jid, "dir": d, "ok": False, "notas": [], "origen": origen}
        tie = clave_valor(os.path.join(d, "tiempo.txt"))
        ent = clave_valor(os.path.join(d, "entorno.txt"))
        c["nodo"] = (ent.get("nodo") or "-").split(".")[0]
        c["cpu"] = ent.get("cpu")
        c["git_hash"] = ent.get("git_hash")
        if not os.path.isdir(d) or not tie.get("estado"):
            st = estado_slurm(jid)
            c["estado"] = st or ("EN_COLA" if not os.path.isdir(d) else "CORRIENDO")
            res.append(c)
            continue
        c.update(leer_outcar(os.path.join(d, "OUTCAR")))
        c["estado"] = tie.get("estado")
        c["wall_s"] = num(tie.get("wall_s"))
        c["afinidad_ok"] = tie.get("afinidad_ok") == "1"
        c["socket"] = tie.get("nucleos_por_socket", "-")
        c["alloc"] = num(tie.get("alloc_cpus"), int)
        c["muestras"] = num(tie.get("hilos_muestras"), int) or 0
        c["nuc_hilos"] = num(tie.get("nucleos_fisicos_hilos"), int)
        c["activos"] = num(tie.get("hilos_activos_mediana"))
        c["uso_cpu"] = num(tie.get("uso_cpu_mediana"))
        c["paso"] = tie.get("paso_variante")
        aj = tie.get("carga_ajena_max", "")
        m = re.match(r"(\d+)/(\d+)", aj)
        c["carga_ajena"] = aj or "-"
        c["carga_frac"] = (int(m.group(1)) / float(m.group(2))) if m and int(m.group(2)) > 0 else None
        if nelm_fijo:
            c["ok"] = (c["estado"] == "OK" and c["afinidad_ok"] and c["wall_s"] is not None
                       and c.get("n_iter") == nelm_fijo)
            if c["estado"] == "OK" and c.get("n_iter") != nelm_fijo:
                c["notas"].append("%s iteraciones (se esperaban %d)" % (c.get("n_iter"), nelm_fijo))
        else:
            c["ok"] = (c["estado"] == "OK" and c["afinidad_ok"] and c.get("convergido") and c["wall_s"] is not None)
        if c["ok"] and c.get("ranks") is not None and c["ranks"] != 1 and origen == "A.2":
            c["ok"] = False
            c["notas"].append("VASP uso %s procesos" % c["ranks"])
        if c["ok"] and origen == "A.2" and c.get("hilos_vasp") is not None and c["hilos_vasp"] != P:
            c["ok"] = False
            c["notas"].append("VASP uso %s hilos" % c["hilos_vasp"])
        if c["ok"] and origen == "A.2" and P > 1 and c["muestras"] and (c["nuc_hilos"] or 0) < P:
            c["ok"] = False
            c["notas"].append("hilos en %s nucleos fisicos" % c["nuc_hilos"])
        if origen == "A.2" and P > 1 and not c["muestras"]:
            c["notas"].append("sin muestras de hilos")
        if jid in exc:
            c["excluida"] = exc[jid]
            c["notas"].append("excluida: " + exc[jid])
        res.append(c)
    return res


# ----------------------------------------------------------------------------- analisis
def analizar(runs):
    avisos = []
    buenas = [c for c in runs if c["ok"] and not c.get("excluida")]
    for campo, nombre in (("n_iter", "iteraciones SCF"), ("nkpts", "NKPTS"), ("nbands", "NBANDS")):
        vals = sorted(set(c.get(campo) for c in buenas if c.get(campo) is not None))
        if len(vals) > 1:
            avisos.append("%s distinto entre corridas: %s (el trabajo no es identico)" % (nombre, ", ".join(str(v) for v in vals)))
    pasos = sorted(set(c.get("paso") for c in buenas if c.get("paso")))
    if len(pasos) > 1:
        avisos.append("las corridas usaron formas distintas de lanzar el paso de srun: %s (ver geometria_pasos.txt)" % ", ".join(pasos))
    cpus = sorted(set(c.get("cpu") for c in buenas if c.get("cpu")))
    if len(cpus) > 1:
        avisos.append("las corridas cayeron en nodos con procesadores distintos: " + " | ".join(cpus))
    es = [float(c["e_sigma0"]) for c in buenas if c.get("e_sigma0")]
    if es and max(es) - min(es) > 1e-3:
        avisos.append("la energia tras 4 iteraciones difiere mas de 1 meV entre corridas (%.6f eV); revisar" % (max(es) - min(es)))
    for P in P_LISTA:
        ts = [c["wall_s"] for c in buenas if c["p"] == P]
        med = mediana(ts)
        for c in buenas:
            if c["p"] != P:
                continue
            if c.get("carga_frac") is not None and c["carga_frac"] >= CARGA_ALTA:
                c["notas"].append("carga ajena alta")
            if med and len(ts) >= 3 and abs(c["wall_s"] - med) / med > DESVIO:
                c["notas"].append("atipica (%+.0f %% de la mediana)" % (100.0 * (c["wall_s"] - med) / med))
    t1, s1 = media_desv([c["wall_s"] for c in buenas if c["p"] == 1])
    i1, si1 = media_desv([c["t_iter_s"] for c in buenas if c["p"] == 1 and c.get("t_iter_s")])
    filas = []
    for P in P_LISTA:
        sel = [c for c in buenas if c["p"] == P]
        tm, ts = media_desv([c["wall_s"] for c in sel])
        it, its = media_desv([c["t_iter_s"] for c in sel if c.get("t_iter_s")])
        mem, _ = media_desv([c["mem_mb"] for c in sel if c.get("mem_mb")])
        f = {"p": P, "n": len(sel), "t_media_s": tm, "t_desv_s": ts,
             "t_min_s": min([c["wall_s"] for c in sel]) if sel else None,
             "t_max_s": max([c["wall_s"] for c in sel]) if sel else None,
             "t_iter_s": it, "t_iter_desv_s": its, "mem_mb": mem,
             "nodos": ",".join(sorted(set(c["nodo"] for c in sel))),
             "sockets": ";".join(sorted(set(c.get("socket", "-") for c in sel)))}
        f["S"], f["S_desv"] = cociente(t1, s1, tm, ts)
        f["S_iter"], f["S_iter_desv"] = cociente(i1, si1, it, its)
        f["E"] = f["S"] / P if f["S"] is not None else None
        f["E_desv"] = f["S_desv"] / P if f["S_desv"] is not None else None
        f["E_iter"] = f["S_iter"] / P if f["S_iter"] is not None else None
        f["E_iter_desv"] = f["S_iter_desv"] / P if f["S_iter_desv"] is not None else None
        filas.append(f)
    return filas, avisos


def comparar_mpi(runs1, filas):
    """S y E por iteracion: MPI (primeras 4 iteraciones de las corridas de la A.1) frente a OpenMP (A.2)"""
    buenas = [c for c in runs1 if c["ok"] and not c.get("excluida") and len(c.get("loops") or []) >= NELM_FIJO]
    por_p = {}
    for c in buenas:
        por_p.setdefault(c["p"], []).append(sum(c["loops"][:NELM_FIJO]) / float(NELM_FIJO))
    m1, s1 = media_desv(por_p.get(1, []))
    comp = []
    for f in filas:
        P = f["p"]
        im, isd = media_desv(por_p.get(P, []))
        S, Sd = cociente(m1, s1, im, isd)
        comp.append({"p": P, "n_mpi": len(por_p.get(P, [])), "t_iter_mpi_s": im, "t_iter_mpi_desv_s": isd,
                     "S_iter_mpi": S, "S_iter_mpi_desv": Sd, "E_iter_mpi": S / P if S else None,
                     "t_iter_omp_s": f["t_iter_s"], "S_iter_omp": f["S_iter"], "S_iter_omp_desv": f["S_iter_desv"],
                     "E_iter_omp": f["E_iter"]})
    return comp


def imprimir(runs, filas, comp, avisos, faltan):
    print("Tarea A.2: escalamiento fuerte con OpenMP")
    print("  " + DESCRIPCION)
    print()
    cab = "%-3s %-3s %-3s %-8s %-10s %-11s %8s %5s %9s %-6s %-7s %-6s %s" % (
        "#", "r", "p", "Job", "Nodo", "Sockets", "t [s]", "Iter", "t/iter", "Hilos", "Activos", "Ajena", "Estado / notas")
    print(cab)
    print("-" * len(cab))
    for c in runs:
        nota = c["estado"] if not c["ok"] else "OK"
        if c["notas"]:
            nota += "; " + "; ".join(c["notas"])
        print("%-3d %-3s %-3d %-8s %-10s %-11s %8s %5s %9s %-6s %-7s %-6s %s" % (
            c["orden"], "r%d" % c["r"], c["p"], c["jobid"], c.get("nodo", "-").replace("granado-", ""),
            c.get("socket", "-"), fmt(c.get("wall_s"), "%.1f"), fmt(c.get("n_iter"), "%d"),
            fmt(c.get("t_iter_s"), "%.1f s"),
            "-" if not c.get("muestras") else "%s/%d" % (c.get("nuc_hilos"), c["p"]),
            "-" if (c.get("activos") is None or not c.get("muestras")) else "%g/%d" % (c["activos"], c["p"]),
            c.get("carga_ajena", "-"), nota))
    print()
    print("t = tiempo de pared de la corrida (arranque + 4 iteraciones + fuerzas y salida). Hilos = nucleos")
    print("fisicos donde corrieron los hilos de VASP. Activos = hilos que usan >= 50 % de un CPU (mediana).")
    print("Ajena = CPUs ocupados por OTROS en los sockets del job.")
    print()
    cab = "%-3s %-2s %15s %13s %13s %13s %13s %13s" % ("p", "n", "t [s]", "S(p)", "E(p)", "t/iter [s]", "S_iter", "E_iter")
    print(cab)
    print("-" * len(cab))
    for f in filas:
        print("%-3d %-2d %15s %13s %13s %13s %13s %13s" % (
            f["p"], f["n"],
            "-" if f["t_media_s"] is None else "%.1f +/- %.1f" % (f["t_media_s"], f["t_desv_s"]),
            "-" if f["S"] is None else "%.2f +/- %.2f" % (f["S"], f["S_desv"]),
            "-" if f["E"] is None else "%.3f +/- %.3f" % (f["E"], f["E_desv"]),
            "-" if f["t_iter_s"] is None else "%.1f +/- %.1f" % (f["t_iter_s"], f["t_iter_desv_s"]),
            "-" if f["S_iter"] is None else "%.2f +/- %.2f" % (f["S_iter"], f["S_iter_desv"]),
            "-" if f["E_iter"] is None else "%.3f +/- %.3f" % (f["E_iter"], f["E_iter_desv"])))
    print()
    print("Por iteracion, MPI (A.1, primeras 4 iteraciones) frente a OpenMP (este barrido):")
    cab = "%-3s %14s %14s %10s %14s %14s %10s" % ("p", "t/iter MPI", "S_iter MPI", "E MPI", "t/iter OMP", "S_iter OMP", "E OMP")
    print(cab)
    print("-" * len(cab))
    for k in comp:
        print("%-3d %14s %14s %10s %14s %14s %10s" % (
            k["p"], fmt(k["t_iter_mpi_s"], "%.1f s"),
            "-" if k["S_iter_mpi"] is None else "%.2f +/- %.2f" % (k["S_iter_mpi"], k["S_iter_mpi_desv"]),
            fmt(k["E_iter_mpi"], "%.3f"), fmt(k["t_iter_omp_s"], "%.1f s"),
            "-" if k["S_iter_omp"] is None else "%.2f +/- %.2f" % (k["S_iter_omp"], k["S_iter_omp_desv"]),
            fmt(k["E_iter_omp"], "%.3f")))
    print()
    print("+/- = desviacion estandar de las repeticiones; la de S y E se propaga de t(1) y t(p).")
    for a in avisos:
        print("AVISO: " + a)
    pend = [c for c in runs if c["estado"] in ("EN_COLA", "CORRIENDO", "PENDING", "RUNNING")]
    print()
    if pend:
        print("Todavia hay %d corridas en cola o corriendo." % len(pend))
    elif faltan:
        print("Faltan corridas validas: %s" % ", ".join(faltan))
        print("Para repetirlas:  bash slurm/parteA2_enviar.sh --faltantes")
    else:
        print("Las corridas de la A.2 estan completas.")


def guardar_csv(runs, filas, comp):
    os.path.isdir(os.path.join(REPO, "datos")) or os.makedirs(os.path.join(REPO, "datos"))
    c1 = ["orden", "p", "r", "jobid", "nodo", "cpu", "estado", "excluida", "socket", "alloc", "afinidad_ok", "paso",
          "muestras", "nuc_hilos", "activos", "uso_cpu", "carga_ajena", "wall_s", "elapsed_s", "n_iter", "t_iter_s",
          "cpu_real", "ranks", "hilos_vasp", "nkpts", "nbands", "mem_mb", "toten", "e_sigma0", "git_hash", "version"]
    with io.open(os.path.join(REPO, "datos", "parteA2_corridas.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(c1)
        for c in runs:
            w.writerow(["" if c.get(k) is None else (("%.4f" % c[k]) if isinstance(c[k], float) else c[k]) for k in c1])
    c2 = ["p", "n", "t_media_s", "t_desv_s", "t_min_s", "t_max_s", "S", "S_desv", "E", "E_desv",
          "t_iter_s", "t_iter_desv_s", "S_iter", "S_iter_desv", "E_iter", "E_iter_desv", "mem_mb", "nodos", "sockets"]
    with io.open(os.path.join(REPO, "datos", "parteA2_resumen.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(c2)
        for f in filas:
            w.writerow(["" if f.get(k) is None else (("%.4f" % f[k]) if isinstance(f[k], float) else f[k]) for k in c2])
    c3 = ["p", "n_mpi", "t_iter_mpi_s", "t_iter_mpi_desv_s", "S_iter_mpi", "S_iter_mpi_desv", "E_iter_mpi",
          "t_iter_omp_s", "S_iter_omp", "S_iter_omp_desv", "E_iter_omp"]
    with io.open(os.path.join(REPO, "datos", "parteA2_mpi_vs_omp.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(c3)
        for k in comp:
            w.writerow(["" if k.get(x) is None else (("%.4f" % k[x]) if isinstance(k[x], float) else k[x]) for x in c3])
    md = [u"# Tarea A.2 - escalamiento fuerte con OpenMP", u"",
          u"Caso base: WSe2 bulk 2x2x1 con V_Se (23 atomos), ENCUT = 400 eV, Gamma 10x10x3 (36 k-points), "
          u"NBANDS = 144. " + DESCRIPCION + u". t = tiempo de pared de `srun` (arranque + 4 iteraciones + fuerzas "
          u"y salida), media +/- desviacion estandar de 3 repeticiones; t(1) con la misma carga.", u"",
          u"| p | n | t (s) | S(p) | E(p) | t/iteracion (s) | S por iteracion | E por iteracion | memoria (MB) |",
          u"|---|---|---|---|---|---|---|---|---|"]
    for f in filas:
        if f["t_media_s"] is None:
            md.append(u"| %d | 0 | - | - | - | - | - | - | - |" % f["p"])
            continue
        md.append(u"| %d | %d | %.1f +/- %.1f | %s | %s | %s | %s | %s | %s |" % (
            f["p"], f["n"], f["t_media_s"], f["t_desv_s"],
            "-" if f["S"] is None else "%.2f +/- %.2f" % (f["S"], f["S_desv"]),
            "-" if f["E"] is None else "%.3f +/- %.3f" % (f["E"], f["E_desv"]),
            "-" if f["t_iter_s"] is None else "%.1f +/- %.1f" % (f["t_iter_s"], f["t_iter_desv_s"]),
            "-" if f["S_iter"] is None else "%.2f +/- %.2f" % (f["S_iter"], f["S_iter_desv"]),
            "-" if f["E_iter"] is None else "%.3f" % f["E_iter"], fmt(f["mem_mb"], "%.0f")))
    md += [u"", u"Por iteracion SCF, MPI (A.1: p procesos con 1 hilo, primeras 4 iteraciones) frente a OpenMP "
                u"(A.2: 1 proceso con p hilos):", u"",
           u"| p | t/iter MPI (s) | S_iter MPI | E_iter MPI | t/iter OpenMP (s) | S_iter OpenMP | E_iter OpenMP |",
           u"|---|---|---|---|---|---|---|"]
    for k in comp:
        md.append(u"| %d | %s | %s | %s | %s | %s | %s |" % (
            k["p"], fmt(k["t_iter_mpi_s"], "%.1f"),
            "-" if k["S_iter_mpi"] is None else "%.2f +/- %.2f" % (k["S_iter_mpi"], k["S_iter_mpi_desv"]),
            fmt(k["E_iter_mpi"], "%.3f"), fmt(k["t_iter_omp_s"], "%.1f"),
            "-" if k["S_iter_omp"] is None else "%.2f +/- %.2f" % (k["S_iter_omp"], k["S_iter_omp_desv"]),
            fmt(k["E_iter_omp"], "%.3f")))
    md.append(u"")
    with io.open(os.path.join(REPO, "datos", "parteA2_tabla.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md))


def guardar_evidencia(runs):
    base = os.path.join(REPO, "datos", "parteA2")
    for c in runs:
        if not os.path.isdir(c["dir"]):
            continue
        dst = os.path.join(base, os.path.basename(c["dir"]))
        os.path.isdir(dst) or os.makedirs(dst)
        for a in ("tiempo.txt", "entorno.txt", "afinidad.txt", "hilos.txt", "hilos_crudo.txt", "progreso.txt", "carga.txt",
                  "OSZICAR", "modulos.txt", "entradas.md5", "slurm_job.txt", "variables.txt", "geometria_pasos.txt"):
            if os.path.isfile(os.path.join(c["dir"], a)):
                shutil.copy(os.path.join(c["dir"], a), os.path.join(dst, a))
        lin = leer(os.path.join(c["dir"], "OUTCAR")).splitlines()
        if lin:
            with io.open(os.path.join(dst, "OUTCAR_extracto.txt"), "w", encoding="utf-8") as fh:
                fh.write(u"\n".join(lin[:300]) + u"\n\n[... recorte ...]\n\n" + u"\n".join(lin[-60:]) + u"\n")
    out = git("add", "datos")
    if out.startswith("ERROR"):
        print("AVISO: git add fallo: " + out)
        return
    if git("diff", "--cached", "--quiet").startswith("ERROR"):
        out = git("commit", "-q", "-m", "Parte A.2: tiempos de escalamiento fuerte OpenMP y evidencia de las corridas")
        print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
    else:
        print("Nada nuevo que guardar en Git.")


def main():
    ap = argparse.ArgumentParser(description="Resumen de la Tarea A.2")
    ap.add_argument("--guardar", action="store_true", help="copia la evidencia a datos/parteA2/ y hace commit")
    a = ap.parse_args()
    runs = corridas(os.path.join(REPO, "logs", "jobs_A2.txt"), "ex_A2",
                    os.path.join(REPO, "datos", "parteA2_excluidas.txt"), "A.2", nelm_fijo=NELM_FIJO)
    if not runs:
        print("Todavia no hay corridas de la A.2: primero  bash slurm/parteA2_enviar.sh")
        return
    runs1 = corridas(os.path.join(REPO, "logs", "jobs_A1.txt"), "ex_A1",
                     os.path.join(REPO, "datos", "parteA1_excluidas.txt"), "A.1")
    filas, avisos = analizar(runs)
    comp = comparar_mpi(runs1, filas)
    if not any(c["ok"] for c in runs1):
        avisos.append("no se encontraron las corridas de la A.1: la comparacion por iteracion con MPI queda vacia")
    faltan = ["p=%d r%d" % (P, R) for P in P_LISTA for R in (1, 2, 3)
              if not any(c["p"] == P and c["r"] == R and c["ok"] and not c.get("excluida") for c in runs)]
    imprimir(runs, filas, comp, avisos, faltan)
    guardar_csv(runs, filas, comp)
    print()
    print("Tablas guardadas en datos/parteA2_corridas.csv, parteA2_resumen.csv, parteA2_mpi_vs_omp.csv y parteA2_tabla.md")
    if a.guardar:
        guardar_evidencia(runs)


if __name__ == "__main__":
    main()
