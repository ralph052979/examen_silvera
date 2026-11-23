#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteB_resumen.py - Parte B: escalamiento debil (B.1) e hibrido MPI x OpenMP (B.2).

B.1  Lee logs/jobs_B1.txt: p procesos MPI (p = 1..24), KPAR = p, 3 k-points
     irreducibles por proceso, SCF completo. Verifica cada corrida (afinidad,
     convergencia, NKPTS = 3p, KPAR = p, p procesos) y reporta t(p), la eficiencia
     debil Ew(p) = t(1)/t(p) (de la corrida y por iteracion SCF) y la memoria
     (VmHWM de cada proceso, medida cada 3 s).
B.2  Lee logs/jobs_B2.txt: 24 nucleos con N procesos MPI x T hilos OpenMP y la carga
     fija de 4 iteraciones SCF del caso base. La combinacion (1,24) suma las corridas
     p = 24 de la A.2 (misma carga, mismo entorno) y t(1) es la p = 1 de la A.2.
     Reporta t, S = t(1)/t, E = S/24, la memoria maxima por proceso y total, los
     procesos que quedaron repartidos en los dos sockets, y la combinacion optima.

Uso (desde la carpeta del repositorio):
  python3 analisis/parteB_resumen.py             # tablas, verificaciones y CSV
  python3 analisis/parteB_resumen.py --guardar   # ademas copia la evidencia a datos/ y hace commit

Archivos que escribe:
  datos/parteB1_corridas.csv   datos/parteB1_resumen.csv
  datos/parteB2_corridas.csv   datos/parteB2_resumen.csv
  datos/parteB_tabla.md        tablas para el reporte

Una corrida se puede excluir agregando su numero de job y el motivo a
datos/parteB_excluidas.txt; luego  bash slurm/parteB_enviar.sh --faltantes  la repite.

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
P_B1 = [1, 2, 4, 8, 12, 16, 24]
MALLAS = {1: "3x3x1", 2: "3x3x2", 4: "5x5x3", 8: "6x6x4", 12: "10x10x3", 16: "14x14x2", 24: "10x10x6"}
COMBOS = [(24, 1), (12, 2), (8, 3), (6, 4), (4, 6), (2, 12), (1, 24)]
NUCLEOS_B2 = 24
NELM_FIJO = 4
DESVIO = 0.15
CARGA_ALTA = 0.25


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
    xs = [x for x in xs if x is not None]
    n = len(xs)
    if n == 0:
        return None, None
    m = sum(xs) / float(n)
    if n == 1:
        return m, 0.0
    return m, math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))


def mediana(xs):
    s = sorted(x for x in xs if x is not None)
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


def excluidas(*rutas):
    d = {}
    for ruta in rutas:
        for l in leer(ruta).splitlines():
            p = l.split(None, 1)
            if p and not p[0].startswith("#"):
                d[p[0]] = p[1].strip() if len(p) > 1 else "sin motivo"
    return d


# ----------------------------------------------------------------------------- lectura
def leer_outcar(ruta):
    t = leer(ruta)
    o = {}
    if not t:
        return o
    m = re.search(r"running\s+(\d+)\s+mpi-ranks,\s+with\s+(\d+)\s+threads/rank", t)
    o["ranks"] = int(m.group(1)) if m else None
    o["hilos_vasp"] = int(m.group(2)) if m else None
    o["kpar"] = buscar(r"distrk:\s+each k-point on\s+\d+\s+cores,\s+(\d+)\s+groups", t, tipo=int)
    o["npar"] = buscar(r"one band on NCORE=\s*\d+\s+cores,\s+(\d+)\s+groups", t, tipo=int)
    o["nkpts"] = buscar(r"NKPTS\s*=\s*(\d+)", t, tipo=int)
    o["nbands"] = buscar(r"number of bands\s+NBANDS\s*=\s*(\d+)", t, ultimo=True, tipo=int)
    pares = [(float(a), float(b)) for a, b in re.findall(r"LOOP:\s+cpu time\s+([\d.]+)\s*:\s*real time\s+([\d.]+)", t)]
    o["loops"] = [b for a, b in pares]
    o["n_iter"] = len(pares)
    o["t_scf_s"] = sum(b for a, b in pares) if pares else None
    o["t_iter_s"] = sum(b for a, b in pares) / len(pares) if pares else None
    o["elapsed_s"] = buscar(r"Elapsed time \(sec\):\s*([\d.]+)", t, tipo=float)
    v = buscar(r"Maximum memory used \(kb\):\s*([\d.]+)", t, tipo=float)
    o["mem_rank0_mb"] = v / 1024.0 if v is not None else None
    o["convergido"] = "aborting loop because EDIFF is reached" in t
    o["toten"] = buscar(r"free\s+energy\s+TOTEN\s*=\s*([-\d.]+)\s*eV", t, ultimo=True)
    o["e_sigma0"] = buscar(r"energy\(sigma->0\)\s*=\s*([-\d.]+)", t, ultimo=True)
    return o


def leer_corrida(c):
    d = c["dir"]
    tie = clave_valor(os.path.join(d, "tiempo.txt"))
    ent = clave_valor(os.path.join(d, "entorno.txt"))
    c["nodo"] = (ent.get("nodo") or "-").split(".")[0]
    c["cpu"] = ent.get("cpu")
    c["git_hash"] = ent.get("git_hash")
    if not os.path.isdir(d) or not tie.get("estado"):
        st = estado_slurm(c["jobid"])
        c["estado"] = st or ("EN_COLA" if not os.path.isdir(d) else "CORRIENDO")
        return False
    c.update(leer_outcar(os.path.join(d, "OUTCAR")))
    c["estado"] = tie.get("estado")
    c["wall_s"] = num(tie.get("wall_s"))
    c["t_fuera_s"] = (c["wall_s"] - c["t_scf_s"]) if (c["wall_s"] is not None and c.get("t_scf_s") is not None) else None
    c["afinidad_ok"] = tie.get("afinidad_ok") == "1"
    c["socket"] = tie.get("nucleos_por_socket", "-")
    c["en_2_sockets"] = num(tie.get("procesos_en_2_sockets"), int)
    c["alloc"] = num(tie.get("alloc_cpus"), int)
    c["muestras"] = num(tie.get("hilos_muestras"), int) or 0
    c["nuc_hilos"] = num(tie.get("nucleos_fisicos_hilos"), int)
    c["activos"] = num(tie.get("hilos_activos_mediana"))
    c["uso_cpu"] = num(tie.get("uso_cpu_mediana"))
    c["paso"] = tie.get("paso_variante")
    c["mem_medidos"] = num(tie.get("mem_procesos_medidos"), int) or 0
    c["mem_proc_mb"] = num(tie.get("mem_proceso_max_mb")) if c["mem_medidos"] else None
    c["mem_total_mb"] = num(tie.get("mem_total_mb")) if c["mem_medidos"] else None
    c["mem_fuente"] = "VmHWM" if c["mem_medidos"] else "-"
    aj = tie.get("carga_ajena_max", "")
    m = re.match(r"(\d+)/(\d+)", aj)
    c["carga_ajena"] = aj or "-"
    c["carga_frac"] = (int(m.group(1)) / float(m.group(2))) if m and int(m.group(2)) > 0 else None
    return True


def registro(ruta, ncampos):
    """lineas del registro con ncampos numeros de configuracion antes del jobid"""
    filas = []
    for orden, l in enumerate(leer(ruta).splitlines(), 1):
        p = l.split()
        if len(p) < ncampos + 1 or not all(x.isdigit() for x in p[:ncampos + 1]):
            continue
        filas.append((orden, [int(x) for x in p[:ncampos]], p[ncampos]))
    return filas


def corridas_b1(exc):
    res = []
    for orden, (P, R), jid in registro(os.path.join(REPO, "logs", "jobs_B1.txt"), 2):
        c = {"parte": "B.1", "origen": "B.1", "orden": orden, "p": P, "r": R, "jobid": jid, "ok": False, "notas": [],
             "dir": os.path.join(REPO, "corridas", "ex_B1_p%d_r%d_%s" % (P, R, jid))}
        if leer_corrida(c):
            c["ok"] = (c["estado"] == "OK" and c["afinidad_ok"] and c.get("convergido") and c["wall_s"] is not None)
            for campo, esperado, nombre in (("ranks", P, "procesos MPI"), ("nkpts", 3 * P, "NKPTS"),
                                            ("kpar", P, "KPAR (grupos de k-points)"), ("hilos_vasp", 1, "hilos por proceso")):
                if c["ok"] and c.get(campo) is not None and c[campo] != esperado:
                    c["ok"] = False
                    c["notas"].append("%s = %s (se esperaba %d)" % (nombre, c[campo], esperado))
            if c["ok"] and c.get("nkpts") is None:
                c["notas"].append("sin NKPTS en el OUTCAR")
        if jid in exc:
            c["excluida"] = exc[jid]
            c["notas"].append("excluida: " + exc[jid])
        res.append(c)
    return res


def verificar_hibrida(c, N, T):
    c["ok"] = (c["estado"] == "OK" and c["afinidad_ok"] and c["wall_s"] is not None and c.get("n_iter") == NELM_FIJO)
    if c["estado"] == "OK" and c.get("n_iter") != NELM_FIJO:
        c["notas"].append("%s iteraciones (se esperaban %d)" % (c.get("n_iter"), NELM_FIJO))
    for campo, esperado, nombre in (("ranks", N, "procesos MPI"), ("hilos_vasp", T, "hilos por proceso")):
        if c["ok"] and c.get(campo) is not None and c[campo] != esperado:
            c["ok"] = False
            c["notas"].append("VASP uso %s = %s (se esperaba %d)" % (nombre, c[campo], esperado))
    if c["ok"] and T > 1 and c["muestras"] and (c["nuc_hilos"] or 0) < N * T:
        c["ok"] = False
        c["notas"].append("hilos en %s nucleos fisicos" % c["nuc_hilos"])
    if T > 1 and not c["muestras"]:
        c["notas"].append("sin muestras de hilos")


def corridas_b2(exc):
    res = []
    for orden, (N, T, R), jid in registro(os.path.join(REPO, "logs", "jobs_B2.txt"), 3):
        c = {"parte": "B.2", "origen": "B.2", "orden": orden, "n": N, "t": T, "r": R, "jobid": jid, "ok": False,
             "notas": [], "dir": os.path.join(REPO, "corridas", "ex_B2_n%dx%d_r%d_%s" % (N, T, R, jid))}
        if leer_corrida(c):
            verificar_hibrida(c, N, T)
        if jid in exc:
            c["excluida"] = exc[jid]
            c["notas"].append("excluida: " + exc[jid])
        res.append(c)
    return res


def corridas_a2(exc, p_lista):
    """corridas de la A.2 (1 proceso x p hilos, misma carga de 4 iteraciones)"""
    res = []
    for orden, (P, R), jid in registro(os.path.join(REPO, "logs", "jobs_A2.txt"), 2):
        if P not in p_lista:
            continue
        c = {"parte": "B.2", "origen": "A.2", "orden": orden, "n": 1, "t": P, "r": R, "jobid": jid, "ok": False,
             "notas": [], "dir": os.path.join(REPO, "corridas", "ex_A2_p%d_r%d_%s" % (P, R, jid))}
        if leer_corrida(c):
            verificar_hibrida(c, 1, P)
            if c["mem_total_mb"] is None and c.get("mem_rank0_mb"):     # 1 proceso: el OUTCAR da su memoria
                c["mem_proc_mb"] = c["mem_total_mb"] = c["mem_rank0_mb"]
                c["mem_fuente"] = "OUTCAR"
        if jid in exc:
            c["excluida"] = exc[jid]
            c["notas"].append("excluida: " + exc[jid])
        res.append(c)
    return res


# ----------------------------------------------------------------------------- analisis
def marcar_atipicas(buenas, clave):
    grupos = {}
    for c in buenas:
        grupos.setdefault(clave(c), []).append(c)
    for g in grupos.values():
        med = mediana([c["wall_s"] for c in g])
        for c in g:
            if c.get("carga_frac") is not None and c["carga_frac"] >= CARGA_ALTA:
                c["notas"].append("carga ajena alta")
            if med and len(g) >= 3 and abs(c["wall_s"] - med) / med > DESVIO:
                c["notas"].append("atipica (%+.0f %% de la mediana)" % (100.0 * (c["wall_s"] - med) / med))


def analizar_b1(runs):
    avisos = []
    buenas = [c for c in runs if c["ok"] and not c.get("excluida")]
    marcar_atipicas(buenas, lambda c: c["p"])
    cpus = sorted(set(c.get("cpu") for c in buenas if c.get("cpu")))
    if len(cpus) > 1:
        avisos.append("B.1: corridas en nodos con procesadores distintos: " + " | ".join(cpus))
    t1, s1 = media_desv([c["wall_s"] for c in buenas if c["p"] == 1])
    i1, si1 = media_desv([c["t_iter_s"] for c in buenas if c["p"] == 1])
    filas = []
    for P in P_B1:
        sel = [c for c in buenas if c["p"] == P]
        tm, ts = media_desv([c["wall_s"] for c in sel])
        it, its = media_desv([c["t_iter_s"] for c in sel])
        es = [float(c["e_sigma0"]) for c in sel if c.get("e_sigma0")]
        f = {"p": P, "malla": MALLAS[P], "nkpts": 3 * P, "n": len(sel), "t_media_s": tm, "t_desv_s": ts,
             "n_iter": media_desv([c["n_iter"] for c in sel])[0],
             "n_iter_lista": "/".join(str(c["n_iter"]) for c in sel),
             "t_iter_s": it, "t_iter_desv_s": its,
             "t_fuera_s": media_desv([c["t_fuera_s"] for c in sel])[0],
             "mem_proc_mb": media_desv([c["mem_proc_mb"] for c in sel])[0],
             "mem_proc_desv_mb": media_desv([c["mem_proc_mb"] for c in sel])[1],
             "mem_total_mb": media_desv([c["mem_total_mb"] for c in sel])[0],
             "mem_rank0_mb": media_desv([c.get("mem_rank0_mb") for c in sel])[0],
             "h_nucleo": sum(P * c["wall_s"] / 3600.0 for c in sel),
             "e_sigma0_eV": media_desv(es)[0], "e_rango_eV": (max(es) - min(es)) if es else None,
             "nodos": ",".join(sorted(set(c["nodo"] for c in sel))),
             "sockets": ";".join(sorted(set(c.get("socket", "-") for c in sel)))}
        f["Ew"], f["Ew_desv"] = cociente(t1, s1, tm, ts)
        f["Ew_iter"], f["Ew_iter_desv"] = cociente(i1, si1, it, its)
        if es and f["e_rango_eV"] > 1e-5:
            avisos.append("B.1 p = %d: la energia difiere %.2e eV entre repeticiones" % (P, f["e_rango_eV"]))
        filas.append(f)
    return filas, avisos


def analizar_b2(runs, runs_t1):
    avisos = []
    buenas = [c for c in runs if c["ok"] and not c.get("excluida")]
    marcar_atipicas(buenas, lambda c: (c["n"], c["t"]))
    uno = [c for c in runs_t1 if c["ok"] and not c.get("excluida")]
    t1, s1 = media_desv([c["wall_s"] for c in uno])
    i1, si1 = media_desv([c["t_iter_s"] for c in uno])
    if not uno:
        avisos.append("B.2: no encontre las corridas p = 1 de la A.2: S y E quedan vacios")
    pasos = sorted(set(c.get("paso") for c in buenas if c.get("paso")))
    if len(pasos) > 1:
        avisos.append("B.2: el paso de srun se lanzo de formas distintas: %s (ver geometria_pasos.txt)" % ", ".join(pasos))
    filas = []
    for N, T in COMBOS:
        sel = [c for c in buenas if (c["n"], c["t"]) == (N, T)]
        tm, ts = media_desv([c["wall_s"] for c in sel])
        it, its = media_desv([c["t_iter_s"] for c in sel])
        f = {"n_mpi": N, "n_omp": T, "n": len(sel),
             "origen": "+".join("%d %s" % (len([c for c in sel if c["origen"] == o]), o)
                                for o in ("B.2", "A.2") if any(c["origen"] == o for c in sel)) or "-",
             "t_media_s": tm, "t_desv_s": ts,
             "t_min_s": min([c["wall_s"] for c in sel]) if sel else None,
             "t_max_s": max([c["wall_s"] for c in sel]) if sel else None,
             "t_iter_s": it, "t_iter_desv_s": its,
             "t_fuera_s": media_desv([c["t_fuera_s"] for c in sel])[0],
             "mem_proc_mb": media_desv([c["mem_proc_mb"] for c in sel])[0],
             "mem_total_mb": media_desv([c["mem_total_mb"] for c in sel])[0],
             "mem_rank0_mb": media_desv([c.get("mem_rank0_mb") for c in sel])[0],
             "mem_fuente": "/".join(sorted(set(c["mem_fuente"] for c in sel if c["mem_fuente"] != "-"))) or "-",
             "en_2_sockets": "/".join(str(c["en_2_sockets"]) for c in sel if c.get("en_2_sockets") is not None) or "-",
             "sockets": ";".join(sorted(set(c.get("socket", "-") for c in sel))),
             "paso": ",".join(sorted(set(c.get("paso") or "-" for c in sel))),
             "h_nucleo": sum(NUCLEOS_B2 * c["wall_s"] / 3600.0 for c in sel if c["origen"] == "B.2"),
             "nodos": ",".join(sorted(set(c["nodo"] for c in sel)))}
        f["S"], f["S_desv"] = cociente(t1, s1, tm, ts)
        f["E"] = f["S"] / NUCLEOS_B2 if f["S"] is not None else None
        f["E_desv"] = f["S_desv"] / NUCLEOS_B2 if f["S_desv"] is not None else None
        f["S_iter"], f["S_iter_desv"] = cociente(i1, si1, it, its)
        filas.append(f)
    nuevas = [c["wall_s"] for c in buenas if (c["n"], c["t"]) == (1, 24) and c["origen"] == "B.2"]
    viejas = [c["wall_s"] for c in buenas if (c["n"], c["t"]) == (1, 24) and c["origen"] == "A.2"]
    if nuevas and viejas:
        a, b = media_desv(nuevas)[0], media_desv(viejas)[0]
        if abs(a - b) / b > DESVIO:
            avisos.append("B.2 (1,24): la corrida nueva (%.0f s) difiere %+.0f %% de las de la A.2 (%.0f s)" % (a, 100.0 * (a - b) / b, b))
    medidas = [f for f in filas if f["t_media_s"] is not None]
    opt = min(medidas, key=lambda f: f["t_media_s"]) if medidas else None
    return filas, avisos, opt, (t1, s1)


# ----------------------------------------------------------------------------- salida
def imprimir_b1(runs, filas):
    print("Parte B.1: escalamiento debil (p procesos MPI, KPAR = p, 3 k-points irreducibles por proceso, SCF completo)")
    print()
    print("Corridas")
    print("%-4s %-3s %-9s %-11s %-10s %7s %5s %5s %5s %8s %9s  %s" % ("p", "r", "Job", "Estado", "Nodo", "t (s)", "Iter",
                                                                     "NKPTS", "KPAR", "Mem/proc", "Mem total", "Notas"))
    for c in sorted(runs, key=lambda c: (c["p"], c["r"], c["orden"])):
        print("%-4d %-3d %-9s %-11s %-10s %7s %5s %5s %5s %8s %9s  %s" % (
            c["p"], c["r"], c["jobid"], ("OK" if c["ok"] else c.get("estado", "?"))[:11], c.get("nodo", "-")[:10],
            fmt(c.get("wall_s"), "%.0f"), fmt(c.get("n_iter"), "%d"), fmt(c.get("nkpts"), "%d"), fmt(c.get("kpar"), "%d"),
            fmt(c.get("mem_proc_mb"), "%.0f"), fmt(c.get("mem_total_mb"), "%.0f"), "; ".join(c["notas"])))
    print("(Mem en MB: maximo VmHWM de un proceso y suma de los maximos de todos los procesos)")
    print()
    print("Resumen (media +/- desviacion estandar de las repeticiones)")
    print("%-3s %-8s %5s %2s %16s %7s %13s %15s %15s %9s %10s" % ("p", "malla", "NKPTS", "n", "t (s)", "Iter", "t/iter (s)",
                                                                   "Ew = t1/t", "Ew por iter", "Mem/proc", "Mem tot GB"))
    for f in filas:
        print("%-3d %-8s %5d %2d %16s %7s %13s %15s %15s %9s %10s" % (
            f["p"], f["malla"], f["nkpts"], f["n"],
            "-" if f["t_media_s"] is None else "%.1f +/- %.1f" % (f["t_media_s"], f["t_desv_s"]),
            fmt(f["n_iter"], "%.1f"),
            "-" if f["t_iter_s"] is None else "%.2f +/- %.2f" % (f["t_iter_s"], f["t_iter_desv_s"]),
            "-" if f["Ew"] is None else "%.3f +/- %.3f" % (f["Ew"], f["Ew_desv"]),
            "-" if f["Ew_iter"] is None else "%.3f +/- %.3f" % (f["Ew_iter"], f["Ew_iter_desv"]),
            fmt(f["mem_proc_mb"], "%.0f"), fmt(f["mem_total_mb"] / 1024.0 if f["mem_total_mb"] else None, "%.2f")))


def imprimir_b2(runs, filas, opt, t1):
    print("Parte B.2: hibrido con 24 nucleos (N procesos MPI x T hilos OpenMP), 4 iteraciones SCF del caso base")
    print("t(1) = %s s (p = 1 de la A.2, misma carga)" % ("-" if t1[0] is None else "%.1f +/- %.1f" % t1))
    print()
    print("Corridas")
    print("%-7s %-4s %-3s %-9s %-11s %-10s %7s %5s %8s %9s %7s %-6s  %s" % (
        "N x T", "Orig", "r", "Job", "Estado", "Nodo", "t (s)", "Iter", "Mem/proc", "Mem total", "2 sock", "Paso", "Notas"))
    for c in sorted(runs, key=lambda c: (COMBOS.index((c["n"], c["t"])) if (c["n"], c["t"]) in COMBOS else 99, c["origen"], c["r"])):
        print("%-7s %-4s %-3d %-9s %-11s %-10s %7s %5s %8s %9s %7s %-6s  %s" % (
            "%dx%d" % (c["n"], c["t"]), c["origen"], c["r"], c["jobid"], ("OK" if c["ok"] else c.get("estado", "?"))[:11],
            c.get("nodo", "-")[:10], fmt(c.get("wall_s"), "%.0f"), fmt(c.get("n_iter"), "%d"),
            fmt(c.get("mem_proc_mb"), "%.0f"), fmt(c.get("mem_total_mb"), "%.0f"), fmt(c.get("en_2_sockets"), "%d"),
            (c.get("paso") or "-")[:6], "; ".join(c["notas"])))
    print("(Mem en MB; en las corridas de la A.2 la memoria es la del OUTCAR, que con 1 proceso es la total.")
    print(" 2 sock = procesos cuya mascara tiene nucleos de los dos sockets)")
    print()
    print("Resumen (media +/- desviacion estandar)")
    print("%-7s %-11s %16s %13s %14s %7s %9s %10s %-8s" % ("N x T", "n", "t (s)", "t/iter (s)", "S = t1/t", "E", "Mem/proc",
                                                           "Mem tot GB", "2 sock"))
    for f in filas:
        print("%-7s %-11s %16s %13s %14s %7s %9s %10s %-8s" % (
            "%dx%d" % (f["n_mpi"], f["n_omp"]), "%d (%s)" % (f["n"], f["origen"]) if f["n"] else "0",
            "-" if f["t_media_s"] is None else "%.1f +/- %.1f" % (f["t_media_s"], f["t_desv_s"]),
            "-" if f["t_iter_s"] is None else "%.2f +/- %.2f" % (f["t_iter_s"], f["t_iter_desv_s"]),
            "-" if f["S"] is None else "%.2f +/- %.2f" % (f["S"], f["S_desv"]),
            fmt(f["E"], "%.3f"), fmt(f["mem_proc_mb"], "%.0f"),
            fmt(f["mem_total_mb"] / 1024.0 if f["mem_total_mb"] else None, "%.2f"), f["en_2_sockets"]))
    if opt:
        print()
        print("Combinacion optima (menor t medio): N x T = %d x %d, t = %.1f +/- %.1f s" % (
            opt["n_mpi"], opt["n_omp"], opt["t_media_s"], opt["t_desv_s"]))


def escribir_csv(ruta, cols, filas):
    with io.open(ruta, "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for f in filas:
            fila = []
            for k in cols:
                v = f.get(k)
                if isinstance(v, bool):
                    v = int(v)
                elif isinstance(v, float):
                    v = "%.12g" % v
                fila.append("" if v is None else v)
            w.writerow(fila)


def guardar_csv(r1, f1, r2, f2, opt, t1):
    d = os.path.join(REPO, "datos")
    for c in r1 + r2:
        c["notas_txt"] = "; ".join(c["notas"])
    escribir_csv(os.path.join(d, "parteB1_corridas.csv"),
                 ["p", "r", "jobid", "estado", "ok", "nodo", "wall_s", "n_iter", "t_iter_s", "t_fuera_s", "ranks", "kpar",
                  "nkpts", "nbands", "mem_proc_mb", "mem_total_mb", "mem_rank0_mb", "socket", "en_2_sockets", "toten",
                  "e_sigma0", "git_hash", "notas_txt"], r1)
    escribir_csv(os.path.join(d, "parteB1_resumen.csv"),
                 ["p", "malla", "nkpts", "n", "t_media_s", "t_desv_s", "n_iter", "n_iter_lista", "t_iter_s", "t_iter_desv_s",
                  "t_fuera_s", "Ew", "Ew_desv", "Ew_iter", "Ew_iter_desv", "mem_proc_mb", "mem_proc_desv_mb", "mem_total_mb",
                  "mem_rank0_mb", "h_nucleo", "e_sigma0_eV", "e_rango_eV", "nodos", "sockets"], f1)
    escribir_csv(os.path.join(d, "parteB2_corridas.csv"),
                 ["n", "t", "origen", "r", "jobid", "estado", "ok", "nodo", "wall_s", "n_iter", "t_iter_s", "t_fuera_s",
                  "ranks", "hilos_vasp", "npar", "mem_proc_mb", "mem_total_mb", "mem_rank0_mb", "mem_fuente", "socket",
                  "en_2_sockets", "paso", "nuc_hilos", "activos", "uso_cpu", "e_sigma0", "git_hash", "notas_txt"], r2)
    escribir_csv(os.path.join(d, "parteB2_resumen.csv"),
                 ["n_mpi", "n_omp", "n", "origen", "t_media_s", "t_desv_s", "t_min_s", "t_max_s", "t_iter_s", "t_iter_desv_s",
                  "t_fuera_s", "S", "S_desv", "E", "E_desv", "S_iter", "S_iter_desv", "mem_proc_mb", "mem_total_mb",
                  "mem_rank0_mb", "mem_fuente", "en_2_sockets", "sockets", "paso", "h_nucleo", "nodos"], f2)
    md = [u"# Parte B - escalamiento debil e hibrido", u"",
          u"## B.1 Escalamiento debil", u"",
          u"p procesos MPI con KPAR = p (un grupo de k-points por proceso) y una malla Gamma con 3p k-points "
          u"irreducibles: cada proceso calcula 3 k-points con las 144 bandas. ENCUT = 400 eV, SCF completo "
          u"(EDIFF = 1E-6). Ew(p) = t(1)/t(p). Memoria: VmHWM de cada proceso de VASP, leida cada 3 s.", u"",
          u"| p | malla | NKPTS | n | t (s) | iter. SCF | t/iter (s) | Ew | Ew por iter. | mem/proceso (MB) | mem total (GB) |",
          u"|---|---|---|---|---|---|---|---|---|---|---|"]
    for f in f1:
        md.append(u"| %d | %s | %d | %d | %s | %s | %s | %s | %s | %s | %s |" % (
            f["p"], f["malla"], f["nkpts"], f["n"],
            "-" if f["t_media_s"] is None else "%.1f +/- %.1f" % (f["t_media_s"], f["t_desv_s"]),
            f["n_iter_lista"] or "-",
            "-" if f["t_iter_s"] is None else "%.2f +/- %.2f" % (f["t_iter_s"], f["t_iter_desv_s"]),
            "-" if f["Ew"] is None else "%.3f +/- %.3f" % (f["Ew"], f["Ew_desv"]),
            "-" if f["Ew_iter"] is None else "%.3f +/- %.3f" % (f["Ew_iter"], f["Ew_iter_desv"]),
            fmt(f["mem_proc_mb"], "%.0f"), fmt(f["mem_total_mb"] / 1024.0 if f["mem_total_mb"] else None, "%.2f")))
    md += [u"", u"## B.2 Hibrido MPI x OpenMP con 24 nucleos", u"",
           u"Carga fija de 4 iteraciones SCF del caso base (NELM = 4); t(1) = %s s es la p = 1 de la A.2 con la misma "
           u"carga. (1,24) incluye las corridas p = 24 de la A.2. Cada proceso con T nucleos fisicos propios "
           u"(mascara verificada), OMP_PLACES=cores, OMP_PROC_BIND=close, OPENBLAS_NUM_THREADS=1." % (
               "-" if t1[0] is None else "%.1f +/- %.1f" % t1), u"",
           u"| N_MPI x N_OMP | n | t (s) | t/iter (s) | S = t(1)/t | E = S/24 | mem/proceso (MB) | mem total (GB) | procesos en 2 sockets |",
           u"|---|---|---|---|---|---|---|---|---|"]
    for f in f2:
        md.append(u"| %d x %d | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            f["n_mpi"], f["n_omp"], "%d (%s)" % (f["n"], f["origen"]) if f["n"] else "0",
            "-" if f["t_media_s"] is None else "%.1f +/- %.1f" % (f["t_media_s"], f["t_desv_s"]),
            "-" if f["t_iter_s"] is None else "%.2f +/- %.2f" % (f["t_iter_s"], f["t_iter_desv_s"]),
            "-" if f["S"] is None else "%.2f +/- %.2f" % (f["S"], f["S_desv"]),
            fmt(f["E"], "%.3f"), fmt(f["mem_proc_mb"], "%.0f"),
            fmt(f["mem_total_mb"] / 1024.0 if f["mem_total_mb"] else None, "%.2f"), f["en_2_sockets"]))
    if opt:
        md += [u"", u"Combinacion optima (menor t medio): **%d x %d**, t = %.1f +/- %.1f s." % (
            opt["n_mpi"], opt["n_omp"], opt["t_media_s"], opt["t_desv_s"])]
    md.append(u"")
    with io.open(os.path.join(d, "parteB_tabla.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md))


def guardar_evidencia(runs):
    base = os.path.join(REPO, "datos", "parteB")
    for c in runs:
        if c.get("origen") == "A.2" or not os.path.isdir(c["dir"]):
            continue
        dst = os.path.join(base, os.path.basename(c["dir"]))
        os.path.isdir(dst) or os.makedirs(dst)
        for a in ("tiempo.txt", "entorno.txt", "afinidad.txt", "hilos.txt", "progreso.txt", "carga.txt", "OSZICAR",
                  "modulos.txt", "entradas.md5", "slurm_job.txt", "variables.txt", "geometria_pasos.txt", "INCAR", "KPOINTS"):
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
        out = git("commit", "-q", "-m", "Parte B: escalamiento debil e hibrido (tiempos, memoria y evidencia de las corridas)")
        print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
    else:
        print("Nada nuevo que guardar en Git.")


def main():
    ap = argparse.ArgumentParser(description="Resumen de la Parte B")
    ap.add_argument("--guardar", action="store_true", help="copia la evidencia a datos/parteB/ y hace commit")
    a = ap.parse_args()
    exc = excluidas(os.path.join(REPO, "datos", "parteB_excluidas.txt"))
    exc_a2 = excluidas(os.path.join(REPO, "datos", "parteA2_excluidas.txt"))
    r1 = corridas_b1(exc)
    r2 = corridas_b2(exc)
    if not r1 and not r2:
        print("Todavia no hay corridas de la Parte B: primero  bash slurm/parteB_enviar.sh")
        return
    a2 = corridas_a2(exc_a2, [1, 24])
    r2 = r2 + [c for c in a2 if c["t"] == 24]
    f1, av1 = analizar_b1(r1)
    f2, av2, opt, t1 = analizar_b2(r2, [c for c in a2 if c["t"] == 1])
    imprimir_b1(r1, f1)
    print()
    print("=" * 100)
    imprimir_b2(r2, f2, opt, t1)
    faltan = (["B.1 p=%d r%d" % (P, R) for P in P_B1 for R in (1, 2, 3)
               if not any(c["p"] == P and c["r"] == R and c["ok"] and not c.get("excluida") for c in r1)] +
              ["B.2 %dx%d r%d" % (N, T, R) for N, T in COMBOS for R in ((1,) if (N, T) == (1, 24) else (1, 2, 3))
               if not any((c["n"], c["t"]) == (N, T) and c["r"] == R and c["origen"] == "B.2" and c["ok"]
                          and not c.get("excluida") for c in r2)])
    avisos = av1 + av2
    print()
    if avisos:
        print("Avisos:")
        for x in avisos:
            print("  * " + x)
    en_cola = [c for c in r1 + r2 if c.get("estado") in ("EN_COLA", "CORRIENDO", "PENDING", "RUNNING", "CONFIGURING", "COMPLETING")]
    if en_cola:
        print("Quedan %d corridas en cola o corriendo: repite este comando cuando terminen." % len(en_cola))
    elif faltan:
        print("Faltan o fallaron: " + ", ".join(faltan))
        print("Mandame una captura antes de repetirlas (bash slurm/parteB_enviar.sh --faltantes).")
    else:
        print("VERIFICACION: las %d corridas de la Parte B estan completas y verificadas." % len([c for c in r1 + r2 if c["origen"] != "A.2" and not c.get("excluida")]))
    guardar_csv(r1, f1, r2, f2, opt, t1)
    print()
    print("Tablas guardadas en datos/parteB1_corridas.csv, parteB1_resumen.csv, parteB2_corridas.csv, parteB2_resumen.csv y parteB_tabla.md")
    if a.guardar:
        guardar_evidencia(r1 + r2)


if __name__ == "__main__":
    main()
