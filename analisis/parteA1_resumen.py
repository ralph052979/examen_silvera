#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteA1_resumen.py - Tarea A.1: escalamiento fuerte con MPI.

Lee las corridas de logs/jobs_A1.txt (21 = 7 valores de p x 3 repeticiones),
verifica cada una y reporta t(p) con media y desviacion estandar. Agrega
S(p) = t(1)/t(p) y E(p) = S(p)/p con sus incertidumbres, que se grafican en la A.3.

Uso (desde la carpeta del repositorio):
  python3 analisis/parteA1_resumen.py             # tablas, verificaciones y CSV
  python3 analisis/parteA1_resumen.py --guardar   # ademas copia la evidencia a datos/ y hace commit

Archivos que escribe:
  datos/parteA1_corridas.csv   una fila por corrida (tiempos crudos, nodo, afinidad, energia...)
  datos/parteA1_resumen.csv    una fila por p (media, desviacion estandar, S, E)
  datos/parteA1_tabla.md       tabla para el reporte

Una corrida se puede excluir del analisis (por ejemplo, si otro usuario saturo el
socket) agregando su numero de job y el motivo a datos/parteA1_excluidas.txt;
luego  bash slurm/parteA1_enviar.sh --faltantes  la repite.

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
REG = os.path.join(REPO, "logs", "jobs_A1.txt")
EXCL = os.path.join(REPO, "datos", "parteA1_excluidas.txt")
P_LISTA = [1, 2, 4, 8, 12, 16, 24]
CARGA_ALTA = 0.25      # fraccion de CPUs logicos de los sockets del job ocupados por otros
DESVIO = 0.15          # una corrida que se aparta mas de 15 % de la mediana de su p es atipica


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


# ----------------------------------------------------------------------------- lectura
def leer_outcar(ruta):
    t = leer(ruta)
    o = {}
    if not t:
        return o
    o["version"] = buscar(r"(vasp\.[0-9.]+[^\n]*)", t)
    o["ranks"] = buscar(r"running\s+(\d+)\s+mpi-ranks", t, tipo=int)
    if o["ranks"] is None:
        o["ranks"] = buscar(r"running on\s+(\d+)\s+total cores", t, tipo=int)
    o["hilos"] = buscar(r"with\s+(\d+)\s+threads/rank", t, tipo=int) or 1
    o["nkpts"] = buscar(r"NKPTS\s*=\s*(\d+)", t, tipo=int)
    o["nbands"] = buscar(r"number of bands\s+NBANDS\s*=\s*(\d+)", t, ultimo=True, tipo=int)
    o["npar"] = buscar(r"one band on NCORE=\s*\d+\s+cores,\s+(\d+)\s+groups", t, tipo=int)
    loops = [float(x) for x in re.findall(r"LOOP:\s+cpu time\s+[\d.]+\s*:\s*real time\s+([\d.]+)", t)]
    o["n_iter"] = len(loops)
    o["t_iter_s"] = sum(loops) / len(loops) if loops else None
    o["elapsed_s"] = buscar(r"Elapsed time \(sec\):\s*([\d.]+)", t, tipo=float)
    v = buscar(r"Maximum memory used \(kb\):\s*([\d.]+)", t, tipo=float)
    o["mem_rank0_mb"] = v / 1024.0 if v is not None else None
    o["convergido"] = "aborting loop because EDIFF is reached" in t
    # energias como texto, con todos los digitos que imprime VASP (para la Tarea C.1)
    o["toten"] = buscar(r"free\s+energy\s+TOTEN\s*=\s*([-\d.]+)\s*eV", t, ultimo=True)
    o["e_sigma0"] = buscar(r"energy\(sigma->0\)\s*=\s*([-\d.]+)", t, ultimo=True)
    return o


def excluidas():
    d = {}
    for l in leer(EXCL).splitlines():
        p = l.split(None, 1)
        if p and not p[0].startswith("#"):
            d[p[0]] = p[1].strip() if len(p) > 1 else "sin motivo"
    return d


def corridas():
    if not os.path.isfile(REG):
        return []
    exc = excluidas()
    res = []
    for orden, l in enumerate(leer(REG).splitlines(), 1):
        p = l.split()
        if len(p) < 3:
            continue
        P, R, jid = int(p[0]), int(p[1]), p[2]
        d = os.path.join(REPO, "corridas", "ex_A1_p%d_r%d_%s" % (P, R, jid))
        c = {"orden": orden, "p": P, "r": R, "jobid": jid, "dir": d, "ok": False, "notas": []}
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
        aj = tie.get("carga_ajena_max", "")
        m = re.match(r"(\d+)/(\d+)", aj)
        c["carga_ajena"] = aj or "-"
        c["carga_frac"] = (int(m.group(1)) / float(m.group(2))) if m and int(m.group(2)) > 0 else None
        c["ok"] = (c["estado"] == "OK" and c["afinidad_ok"] and c.get("convergido") and c["wall_s"] is not None)
        if c["ok"] and c.get("ranks") is not None and c["ranks"] != P:
            c["ok"] = False
            c["notas"].append("VASP uso %s procesos" % c["ranks"])
        if jid in exc:
            c["excluida"] = exc[jid]
            c["notas"].append("excluida: " + exc[jid])
        res.append(c)
    return res


# ----------------------------------------------------------------------------- analisis
def analizar(runs):
    avisos = []
    buenas = [c for c in runs if c["ok"] and not c.get("excluida")]
    # mismo trabajo en todas: iteraciones SCF, NKPTS y NBANDS
    for campo, nombre in (("n_iter", "iteraciones SCF"), ("nkpts", "NKPTS"), ("nbands", "NBANDS")):
        vals = sorted(set(c.get(campo) for c in buenas if c.get(campo) is not None))
        if len(vals) > 1:
            avisos.append("%s distinto entre corridas: %s (el trabajo no es identico; ver t por iteracion)"
                          % (nombre, ", ".join(str(v) for v in vals)))
    cpus = sorted(set(c.get("cpu") for c in buenas if c.get("cpu")))
    if len(cpus) > 1:
        avisos.append("las corridas cayeron en nodos con procesadores distintos: " + " | ".join(cpus))
    es = [float(c["e_sigma0"]) for c in buenas if c.get("e_sigma0")]
    if es and max(es) - min(es) > 1e-3:
        avisos.append("la energia difiere mas de 1 meV entre corridas (%.6f eV); revisar" % (max(es) - min(es)))
    # carga ajena y atipicas
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
    # resumen por p
    filas = []
    t1, s1 = media_desv([c["wall_s"] for c in buenas if c["p"] == 1])
    for P in P_LISTA:
        sel = [c for c in buenas if c["p"] == P]
        tm, ts = media_desv([c["wall_s"] for c in sel])
        it, _ = media_desv([c["t_iter_s"] for c in sel if c.get("t_iter_s")])
        f = {"p": P, "n": len(sel), "t_media_s": tm, "t_desv_s": ts,
             "t_min_s": min([c["wall_s"] for c in sel]) if sel else None,
             "t_max_s": max([c["wall_s"] for c in sel]) if sel else None,
             "t_iter_s": it, "S": None, "S_desv": None, "E": None, "E_desv": None,
             "nodos": ",".join(sorted(set(c["nodo"] for c in sel))),
             "sockets": ";".join(sorted(set(c.get("socket", "-") for c in sel)))}
        if tm and t1:
            S = t1 / tm
            rel = math.sqrt((s1 / t1) ** 2 + (ts / tm) ** 2)
            f.update({"S": S, "S_desv": S * rel, "E": S / P, "E_desv": S * rel / P})
        filas.append(f)
    faltan = ["p=%d r%d" % (P, R) for P in P_LISTA for R in (1, 2, 3)
              if not any(c["p"] == P and c["r"] == R and c["ok"] and not c.get("excluida") for c in runs)]
    return filas, avisos, faltan


def imprimir(runs, filas, avisos, faltan):
    print("Tarea A.1: escalamiento fuerte con MPI (caso base, NCORE = 1, KPAR = 1)")
    print()
    cab = "%-3s %-3s %-3s %-8s %-10s %-11s %9s %5s %8s %-6s %-6s %s" % (
        "#", "r", "p", "Job", "Nodo", "Sockets", "t [s]", "Iter", "t/iter", "Afin.", "Ajena", "Estado / notas")
    print(cab)
    print("-" * len(cab))
    for c in runs:
        nota = c["estado"] if not c["ok"] else "OK"
        if c["notas"]:
            nota += "; " + "; ".join(c["notas"])
        print("%-3d %-3s %-3d %-8s %-10s %-11s %9s %5s %8s %-6s %-6s %s" % (
            c["orden"], "r%d" % c["r"], c["p"], c["jobid"], c.get("nodo", "-").replace("granado-", ""),
            c.get("socket", "-"),
            fmt(c.get("wall_s"), "%.1f"), fmt(c.get("n_iter"), "%d"), fmt(c.get("t_iter_s"), "%.2f s"),
            ("ok" if c.get("afinidad_ok") else ("MAL" if "afinidad_ok" in c else "-")),
            c.get("carga_ajena", "-"), nota))
    print()
    print("Sockets = nucleos fisicos usados en cada socket (s0, s1). Ajena = maximo de CPUs logicos")
    print("ocupados por OTROS procesos en los sockets del job durante la corrida.")
    print()
    cab = "%-3s %-2s %16s %8s %17s %15s %15s" % ("p", "n", "t medio [s]", "CV [%]", "min - max [s]", "S(p)", "E(p)")
    print(cab)
    print("-" * len(cab))
    for f in filas:
        cv = 100.0 * f["t_desv_s"] / f["t_media_s"] if f["t_media_s"] else None
        print("%-3d %-2d %16s %8s %17s %15s %15s" % (
            f["p"], f["n"],
            "-" if f["t_media_s"] is None else "%.1f +/- %.1f" % (f["t_media_s"], f["t_desv_s"]),
            fmt(cv, "%.1f"),
            "-" if f["t_min_s"] is None else "%.1f - %.1f" % (f["t_min_s"], f["t_max_s"]),
            "-" if f["S"] is None else "%.2f +/- %.2f" % (f["S"], f["S_desv"]),
            "-" if f["E"] is None else "%.3f +/- %.3f" % (f["E"], f["E_desv"])))
    print()
    print("t = tiempo de pared de srun (date antes y despues); +/- = desviacion estandar de las")
    print("repeticiones; la incertidumbre de S y E se propaga de las desviaciones de t(1) y t(p).")
    for a in avisos:
        print("AVISO: " + a)
    pendientes = [c for c in runs if c["estado"] in ("EN_COLA", "CORRIENDO", "PENDING", "RUNNING")]
    if pendientes:
        print()
        print("Todavia hay %d corridas en cola o corriendo." % len(pendientes))
    elif faltan:
        print()
        print("Faltan corridas validas: %s" % ", ".join(faltan))
        print("Para repetirlas:  bash slurm/parteA1_enviar.sh --faltantes")
    else:
        print()
        print("Las 21 corridas estan completas.")


def guardar_csv(runs, filas):
    os.path.isdir(os.path.join(REPO, "datos")) or os.makedirs(os.path.join(REPO, "datos"))
    c1 = ["orden", "p", "r", "jobid", "nodo", "cpu", "estado", "excluida", "socket", "afinidad_ok", "carga_ajena",
          "wall_s", "elapsed_s", "n_iter", "t_iter_s", "ranks", "hilos", "nkpts", "nbands", "npar", "mem_rank0_mb",
          "toten", "e_sigma0", "git_hash", "version"]
    with io.open(os.path.join(REPO, "datos", "parteA1_corridas.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(c1)
        for c in runs:
            w.writerow(["" if c.get(k) is None else (("%.4f" % c[k]) if isinstance(c[k], float) else c[k])
                        for k in c1])
    c2 = ["p", "n", "t_media_s", "t_desv_s", "t_min_s", "t_max_s", "t_iter_s", "S", "S_desv", "E", "E_desv",
          "nodos", "sockets"]
    with io.open(os.path.join(REPO, "datos", "parteA1_resumen.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(c2)
        for f in filas:
            w.writerow(["" if f.get(k) is None else (("%.4f" % f[k]) if isinstance(f[k], float) else f[k]) for k in c2])
    md = [u"# Tarea A.1 - escalamiento fuerte con MPI", u"",
          u"Caso base: WSe2 bulk 2x2x1 con V_Se (23 atomos), ENCUT = 400 eV, Gamma 10x10x3 (36 k-points), "
          u"NBANDS = 144, NCORE = 1, KPAR = 1. Un nucleo fisico por proceso MPI (`--hint=nomultithread`). "
          u"t = tiempo de pared de `srun`, media +/- desviacion estandar de las repeticiones.", u"",
          u"| p | n | t medio (s) | t medio (min) | min - max (s) | S(p) | E(p) |", u"|---|---|---|---|---|---|---|"]
    for f in filas:
        if f["t_media_s"] is None:
            md.append(u"| %d | 0 | - | - | - | - | - |" % f["p"])
            continue
        md.append(u"| %d | %d | %.1f +/- %.1f | %.2f | %.1f - %.1f | %s | %s |" % (
            f["p"], f["n"], f["t_media_s"], f["t_desv_s"], f["t_media_s"] / 60.0, f["t_min_s"], f["t_max_s"],
            "-" if f["S"] is None else "%.2f +/- %.2f" % (f["S"], f["S_desv"]),
            "-" if f["E"] is None else "%.3f +/- %.3f" % (f["E"], f["E_desv"])))
    md.append(u"")
    with io.open(os.path.join(REPO, "datos", "parteA1_tabla.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md))


def guardar_evidencia(runs):
    base = os.path.join(REPO, "datos", "parteA1")
    for c in runs:
        if not os.path.isdir(c["dir"]):
            continue
        dst = os.path.join(base, os.path.basename(c["dir"]))
        os.path.isdir(dst) or os.makedirs(dst)
        for a in ("tiempo.txt", "entorno.txt", "afinidad.txt", "carga.txt", "OSZICAR", "modulos.txt", "entradas.md5"):
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
        out = git("commit", "-q", "-m", "Parte A.1: tiempos de escalamiento fuerte MPI y evidencia de las corridas")
        print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
    else:
        print("Nada nuevo que guardar en Git.")


def main():
    ap = argparse.ArgumentParser(description="Resumen de la Tarea A.1")
    ap.add_argument("--guardar", action="store_true", help="copia la evidencia a datos/parteA1/ y hace commit")
    a = ap.parse_args()
    runs = corridas()
    if not runs:
        print("Todavia no hay corridas de la A.1: primero  bash slurm/parteA1_enviar.sh")
        return
    filas, avisos, faltan = analizar(runs)
    imprimir(runs, filas, avisos, faltan)
    guardar_csv(runs, filas)
    print()
    print("Tablas guardadas en datos/parteA1_corridas.csv, datos/parteA1_resumen.csv y datos/parteA1_tabla.md")
    if a.guardar:
        guardar_evidencia(runs)


if __name__ == "__main__":
    main()
