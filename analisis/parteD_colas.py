#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteD_colas.py - Parte D: consumo de horas-nucleo, tiempos de cola y costo por calculo.

Lee datos/sacct_colas.txt (salida de sacct con -P). Si no existe y sacct esta disponible
(en Granado), lo genera con:
  sacct -u $USER -S 2026-11-16 -X -P --format=JobID,JobName%40,Partition,AllocCPUS,ReqMem,
        Timelimit,Submit,Eligible,Start,End,Elapsed,CPUTimeRAW,State,NodeList
y guarda tambien el registro de consumo que pide el examen (datos/sacct_registro_examen.txt):
  sacct -u $USER -S 2026-09-07 -X --format=JobID,JobName,Partition,AllocCPUS,Elapsed,CPUTimeRAW,MaxRSS,State

Calcula:
  * horas-nucleo por parte del examen (por el nombre del job)
  * espera en cola de cada job = Start - Eligible (Eligible es cuando la dependencia se cumple)
  * costo en horas-nucleo y tiempo hasta el resultado de un calculo del caso base (SCF completo)
    con cada configuracion medida (A.1 y la B.1 con p = 12, que es el caso base con KPAR = 12)

Escribe: datos/parteD_colas.csv, datos/parteD_consumo.csv, datos/parteD_costo.csv, datos/parteD_tabla.md
Uso:     python3 analisis/parteD_colas.py [--guardar]
Todo lo que imprime es ASCII.
"""
from __future__ import print_function, division
import argparse
import csv
import datetime
import io
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(REPO, "datos")
SACCT = os.path.join(DATOS, "sacct_colas.txt")
REGISTRO = os.path.join(DATOS, "sacct_registro_examen.txt")
PRESUPUESTO = 250.0
CATEGORIAS = [
    ("Tarea 0: calibracion final (K1-K7)", r"^ex_T0b_"),
    ("A.1: escalamiento fuerte MPI", r"^ex_A1_"),
    ("A.2: diagnostico de OpenMP (2 rondas)", r"^ex_A2[dh]_"),
    ("A.2: barrido OpenMP (4 iteraciones)", r"^ex_A2_p"),
    ("B.1: escalamiento debil", r"^ex_B1_"),
    ("B.2: hibrido MPI x OpenMP", r"^ex_B2_"),
    ("C.2: serie de ENCUT", r"^ex_C2_"),
    ("Pruebas e intentos previos (17-18 nov)", r".*"),
]


def git(*args):
    try:
        return subprocess.check_output(("git",) + args, cwd=REPO, stderr=subprocess.STDOUT).decode().strip()
    except (OSError, subprocess.CalledProcessError) as e:
        out = getattr(e, "output", b"")
        return "ERROR: " + (out.decode().strip() if out else str(e))


def generar_sacct():
    usuario = os.environ.get("USER") or ""
    try:
        out = subprocess.check_output(
            ["sacct", "-u", usuario, "-S", "2026-11-16", "-X", "-P",
             "--format=JobID,JobName%40,Partition,AllocCPUS,ReqMem,Timelimit,Submit,Eligible,Start,End,Elapsed,"
             "CPUTimeRAW,State,NodeList"], stderr=subprocess.STDOUT).decode()
        with io.open(SACCT, "w", encoding="utf-8") as fh:
            fh.write(out)
        return True
    except (OSError, subprocess.CalledProcessError):
        return False


def generar_registro():
    """registro de consumo con el comando exacto del enunciado (solo donde hay sacct)"""
    usuario = os.environ.get("USER") or ""
    try:
        reg = subprocess.check_output(
            ["sacct", "-u", usuario, "-S", "2026-09-07", "-X",
             "--format=JobID,JobName,Partition,AllocCPUS,Elapsed,CPUTimeRAW,MaxRSS,State"],
            stderr=subprocess.STDOUT).decode()
    except (OSError, subprocess.CalledProcessError):
        return False
    with io.open(REGISTRO, "w", encoding="utf-8") as fh:
        fh.write(u"$ sacct -u $USER -S 2026-09-07 -X --format=JobID,JobName,Partition,AllocCPUS,Elapsed,CPUTimeRAW,MaxRSS,State\n")
        fh.write(reg)
    return True


def fecha(s):
    try:
        return datetime.datetime.strptime(s.strip(), "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        return None


def percentil(xs, q):
    s = sorted(xs)
    if not s:
        return None
    k = (len(s) - 1) * q
    i = int(k)
    j = min(i + 1, len(s) - 1)
    return s[i] + (s[j] - s[i]) * (k - i)


def leer_csv(nombre):
    ruta = os.path.join(DATOS, nombre)
    if not os.path.isfile(ruta):
        return []
    with io.open(ruta, "r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser(description="Consumo, colas y costo por calculo (Parte D)")
    ap.add_argument("--guardar", action="store_true", help="hace commit de las tablas")
    a = ap.parse_args()
    if not os.path.isfile(SACCT) and not generar_sacct():
        print("[FALLA] no encuentro datos/sacct_colas.txt y no pude correr sacct")
        sys.exit(1)
    if generar_registro():
        print("Registro de consumo del enunciado guardado en datos/sacct_registro_examen.txt")
    with io.open(SACCT, "r", encoding="utf-8") as fh:
        filas = [l.rstrip("\n").split("|") for l in fh if l.strip()]
    cab = filas[0]
    idx = dict((c, i) for i, c in enumerate(cab))
    jobs = []
    for f in filas[1:]:
        if len(f) < len(cab):
            continue
        j = dict((c, f[idx[c]]) for c in cab)
        if j["Partition"] not in ("phys_fast", "phys_general"):
            continue
        j["h"] = float(j["CPUTimeRAW"] or 0) / 3600.0
        e, s = fecha(j["Eligible"]), fecha(j["Start"])
        j["espera_s"] = (s - e).total_seconds() if (e and s) else None
        for nombre, patron in CATEGORIAS:
            if re.match(patron, j["JobName"]):
                j["categoria"] = nombre
                break
        jobs.append(j)
    total = sum(j["h"] for j in jobs)
    consumo = []
    for nombre, patron in CATEGORIAS:
        sel = [j for j in jobs if j["categoria"] == nombre]
        consumo.append({"categoria": nombre, "jobs": len(sel), "h_nucleo": sum(j["h"] for j in sel)})
    examen = [j for j in jobs if j["JobName"].startswith("ex_") and j["espera_s"] is not None]
    esperas = [j["espera_s"] for j in examen]
    inmediatos = len([w for w in esperas if w < 60])
    print("Consumo en phys_fast y phys_general desde 2026-11-16: %d jobs, %.2f h-nucleo de %.0f (%.1f %%)" % (
        len(jobs), total, PRESUPUESTO, 100.0 * total / PRESUPUESTO))
    for c in consumo:
        print("  %-42s %4d jobs  %7.2f h-nucleo" % (c["categoria"], c["jobs"], c["h_nucleo"]))
    print()
    print("Espera en cola de los jobs del examen (Start - Eligible): %d jobs" % len(esperas))
    print("  mediana %.0f s, percentil 90 %.0f s, maximo %.0f s; %d de %d (%.0f %%) empezaron en menos de 1 min" % (
        percentil(esperas, 0.5), percentil(esperas, 0.9), max(esperas), inmediatos, len(esperas),
        100.0 * inmediatos / len(esperas)))
    largos = sorted([j for j in examen if j["espera_s"] >= 60], key=lambda j: -j["espera_s"])[:8]
    for j in largos:
        print("  espera %6.1f min: %s (%s CPUs, elegible %s)" % (j["espera_s"] / 60.0, j["JobName"], j["AllocCPUS"],
                                                                j["Eligible"]))
    por_cpus = {}
    for j in examen:
        k = "1" if j["AllocCPUS"] == "1" else ("2-12" if int(j["AllocCPUS"]) <= 12 else "16-24")
        por_cpus.setdefault(k, []).append(j["espera_s"])
    for k in ("1", "2-12", "16-24"):
        if k in por_cpus:
            print("  %-6s CPUs: %3d jobs, mediana %.0f s, maximo %.0f s" % (k, len(por_cpus[k]), percentil(por_cpus[k], 0.5),
                                                                          max(por_cpus[k])))
    # costo de un calculo del caso base (SCF completo) con cada configuracion medida
    costo = []
    for f in leer_csv("parteA1_resumen.csv"):
        p, t, sd = int(float(f["p"])), float(f["t_media_s"]), float(f["t_desv_s"])
        costo.append({"configuracion": "%d procesos MPI, KPAR = 1 (A.1)" % p, "nucleos": p, "t_s": t, "t_desv_s": sd,
                      "h_nucleo": p * t / 3600.0})
    for f in leer_csv("parteB1_resumen.csv"):
        if f["p"] == "12":
            t, sd = float(f["t_media_s"]), float(f["t_desv_s"])
            costo.append({"configuracion": "12 procesos MPI, KPAR = 12 (B.1, malla del caso base)", "nucleos": 12,
                          "t_s": t, "t_desv_s": sd, "h_nucleo": 12 * t / 3600.0})
    med = percentil(esperas, 0.5) or 0.0
    for c in costo:
        c["t_resultado_min"] = (c["t_s"] + med) / 60.0
    print()
    print("Costo de un calculo del caso base (SCF completo, 18 iteraciones):")
    for c in costo:
        print("  %-55s t = %6.1f +/- %4.1f s  -> %.2f h-nucleo" % (c["configuracion"], c["t_s"], c["t_desv_s"], c["h_nucleo"]))
    with io.open(os.path.join(DATOS, "parteD_colas.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["JobID", "JobName", "categoria", "AllocCPUS", "Submit", "Eligible", "Start", "Elapsed", "h_nucleo",
                    "espera_s", "State"])
        for j in jobs:
            w.writerow([j["JobID"], j["JobName"], j["categoria"], j["AllocCPUS"], j["Submit"], j["Eligible"], j["Start"],
                        j["Elapsed"], "%.4f" % j["h"], "" if j["espera_s"] is None else "%.0f" % j["espera_s"], j["State"]])
    with io.open(os.path.join(DATOS, "parteD_consumo.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["categoria", "jobs", "h_nucleo"])
        for c in consumo:
            w.writerow([c["categoria"], c["jobs"], "%.2f" % c["h_nucleo"]])
        w.writerow(["TOTAL", len(jobs), "%.2f" % total])
    with io.open(os.path.join(DATOS, "parteD_costo.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["configuracion", "nucleos", "t_s", "t_desv_s", "h_nucleo", "t_resultado_min"])
        for c in costo:
            w.writerow([c["configuracion"], c["nucleos"], "%.1f" % c["t_s"], "%.1f" % c["t_desv_s"], "%.3f" % c["h_nucleo"],
                        "%.2f" % c["t_resultado_min"]])
    md = [u"# Parte D - consumo, colas y costo por calculo", u"",
          u"Consumo en phys_fast y phys_general desde 2026-11-16 (sacct, CPUTimeRAW): **%.2f h-nucleo** de %.0f." % (
              total, PRESUPUESTO), u"", u"| parte | jobs | h-nucleo |", u"|---|---|---|"]
    md += [u"| %s | %d | %.2f |" % (c["categoria"], c["jobs"], c["h_nucleo"]) for c in consumo]
    md += [u"| **total** | %d | **%.2f** |" % (len(jobs), total), u"",
           u"Espera en cola de los %d jobs del examen (Start - Eligible): mediana %.0f s, percentil 90 %.0f s, maximo %.1f min; "
           u"%.0f %% empezaron en menos de 1 minuto." % (len(esperas), percentil(esperas, 0.5), percentil(esperas, 0.9),
                                                       max(esperas) / 60.0, 100.0 * inmediatos / len(esperas)), u"",
           u"| configuracion | nucleos | t (s) | h-nucleo por calculo |", u"|---|---|---|---|"]
    md += [u"| %s | %d | %.1f +/- %.1f | %.2f |" % (c["configuracion"], c["nucleos"], c["t_s"], c["t_desv_s"], c["h_nucleo"])
           for c in costo]
    md.append(u"")
    with io.open(os.path.join(DATOS, "parteD_tabla.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md))
    print()
    print("Tablas guardadas en datos/parteD_colas.csv, parteD_consumo.csv, parteD_costo.csv y parteD_tabla.md")
    if a.guardar:
        out = git("add", "datos", "analisis/parteD_colas.py")
        if out.startswith("ERROR"):
            print("AVISO: git add fallo: " + out)
            return
        if git("diff", "--cached", "--quiet").startswith("ERROR"):
            out = git("commit", "-q", "-m", "Parte D: consumo por parte, tiempos de cola y costo por calculo")
            print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
        else:
            print("Nada nuevo que guardar en Git.")


if __name__ == "__main__":
    main()
