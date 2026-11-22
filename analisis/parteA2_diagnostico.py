#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteA2_diagnostico.py - Tarea A.2, paso 1: resultado del diagnostico de OpenMP.

Para cada variante (A, B, C, D) y cada p (4 y 24) dice si los p hilos de VASP
corrieron en p nucleos FISICOS distintos, cuantas CPUs cobro Slurm y si cabe en
la QOS. Elige la variante para el barrido y la escribe en datos/parteA2_variante.txt.

Uso (desde la carpeta del repositorio):
  python3 analisis/parteA2_diagnostico.py             # tabla y recomendacion
  python3 analisis/parteA2_diagnostico.py --guardar   # ademas guarda la evidencia en datos/ y hace commit

Solo usa la biblioteca estandar de Python. Todo lo que imprime es ASCII.
"""
from __future__ import print_function
import argparse
import csv
import io
import os
import re
import shutil
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REG = os.path.join(REPO, "logs", "jobs_A2d.txt")
SALIDA = os.path.join(REPO, "datos", "parteA2_variante.txt")
LIMITE_QOS = 24
DESCRIPCION = {
    "A": "--ntasks=1 --cpus-per-task=p --hint=nomultithread",
    "B": "--ntasks=1 --cpus-per-task=p --threads-per-core=1",
    "C": "--ntasks=p --cpus-per-task=1 --hint=nomultithread + srun -n 1 -c p",
    "D": "--ntasks=1 --cpus-per-task=2p --hint=multithread",
}


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


def num(x, tipo=int):
    try:
        return tipo(x)
    except (TypeError, ValueError):
        return None


def git(*args):
    try:
        return subprocess.check_output(("git",) + args, cwd=REPO, stderr=subprocess.STDOUT).decode().strip()
    except (OSError, subprocess.CalledProcessError) as e:
        out = getattr(e, "output", b"")
        return "ERROR: " + (out.decode().strip() if out else str(e))


def slurm(jobid):
    """(estado, motivo) segun squeue; si ya no esta en cola, el estado de sacct"""
    try:
        out = subprocess.check_output(["squeue", "-h", "-j", str(jobid), "-o", "%T|%r"],
                                      stderr=subprocess.STDOUT).decode().strip()
        if out and "|" in out:
            st, mot = out.splitlines()[0].split("|", 1)
            return st, mot
    except (OSError, subprocess.CalledProcessError):
        pass
    try:
        out = subprocess.check_output(["sacct", "-n", "-X", "-j", str(jobid), "-o", "State%30"],
                                      stderr=subprocess.STDOUT).decode().split()
        if out:
            return out[0], ""
    except (OSError, subprocess.CalledProcessError):
        pass
    return None, ""


def leer_outcar(ruta):
    t = leer(ruta)
    o = {}
    m = re.search(r"running\s+(\d+)\s+mpi-ranks,\s+with\s+(\d+)\s+threads/rank", t)
    o["ranks"] = int(m.group(1)) if m else None
    o["hilos_vasp"] = int(m.group(2)) if m else None
    loops = [float(x) for x in re.findall(r"LOOP:\s+cpu time\s+[\d.]+\s*:\s*real time\s+([\d.]+)", t)]
    o["n_iter"] = len(loops)
    o["t_iter_s"] = sum(loops) / len(loops) if loops else None
    return o


def evaluar():
    filas = []
    for l in leer(REG).splitlines():
        p = l.split()
        if len(p) < 3:
            continue
        V, P, jid = p[0], int(p[1]), p[2]
        f = {"variante": V, "p": P, "jobid": jid, "veredicto": "", "motivo": "", "costo": None}
        if jid == "RECHAZADO":
            f["veredicto"] = "RECHAZADO"
            f["motivo"] = " ".join(p[3:])[:110]
            filas.append(f)
            continue
        d = os.path.join(REPO, "corridas", "ex_A2d_%s%d_%s" % (V, P, jid))
        f["dir"] = d
        tie = clave_valor(os.path.join(d, "tiempo.txt"))
        if not tie.get("estado"):
            st, mot = slurm(jid)
            if st in ("PENDING", "CONFIGURING"):
                f["veredicto"] = "EN COLA"
                f["motivo"] = mot
                if re.search(r"Max(Cpu|TRES)", mot or "", re.I):
                    f["veredicto"] = "NO CABE"
                    f["motivo"] = "Slurm no lo correra nunca (%s): cancelalo con  scancel %s" % (mot, jid)
            elif st in ("RUNNING", "COMPLETING", None):
                f["veredicto"] = "CORRIENDO"
            else:
                f["veredicto"] = "NO PASA"
                f["motivo"] = "el job termino como %s sin resultados" % st
            filas.append(f)
            continue
        o = leer_outcar(os.path.join(d, "OUTCAR"))
        f.update(o)
        f["estado"] = tie.get("estado")
        f["alloc"] = num(tie.get("alloc_cpus"))
        f["nuc_mascara"] = num(tie.get("nucleos_fisicos_distintos"))
        f["nuc_hilos"] = num(tie.get("nucleos_fisicos_hilos"))
        f["muestras"] = num(tie.get("hilos_muestras")) or 0
        f["socket"] = tie.get("nucleos_por_socket", "-")
        if f["alloc"]:
            f["costo"] = f["alloc"] / float(P)
        problemas = []
        if (f["estado"] or "").startswith("FALLA_srun"):
            err = [l.strip() for l in leer(os.path.join(d, "vasp.err")).splitlines() if "error" in l.lower()]
            problemas.append("srun no pudo lanzar VASP" + (": " + err[0][:80] if err else " (%s)" % f["estado"]))
        elif f["nuc_mascara"] is not None and f["nuc_mascara"] < P:
            problemas.append("la mascara tiene %d nucleos fisicos de %d" % (f["nuc_mascara"], P))
        if f["muestras"] and f["nuc_hilos"] is not None and f["nuc_hilos"] < P:
            problemas.append("los hilos corrieron en %d nucleos fisicos de %d" % (f["nuc_hilos"], P))
        if f.get("hilos_vasp") is not None and f["hilos_vasp"] != P:
            problemas.append("VASP uso %d hilos" % f["hilos_vasp"])
        if f.get("ranks") is not None and f["ranks"] != 1:
            problemas.append("VASP uso %d procesos MPI" % f["ranks"])
        if f["estado"] != "OK" and not problemas:
            problemas.append("estado %s" % f["estado"])
        if problemas:
            f["veredicto"] = "NO PASA"
            f["motivo"] = "; ".join(problemas)
        elif f["muestras"] == 0:
            f["veredicto"] = "PASA*"
            f["motivo"] = "mascara correcta, pero sin muestras de los hilos"
        else:
            f["veredicto"] = "PASA"
            f["motivo"] = "%d hilos en %d nucleos fisicos distintos" % (P, f["nuc_hilos"])
        if f["veredicto"].startswith("PASA") and f["costo"] and f["costo"] > 1.01:
            f["motivo"] += "; Slurm cobra %d CPUs (x%.0f)" % (f["alloc"], f["costo"])
        filas.append(f)
    return filas


def recomendar(filas):
    def ver(V, P):
        for f in filas:
            if f["variante"] == V and f["p"] == P:
                return f
        return None

    def pasa(f):
        return f is not None and f["veredicto"].startswith("PASA")

    for V in ("C", "B"):
        f4, f24 = ver(V, 4), ver(V, 24)
        if pasa(f4) and pasa(f24) and (f24["costo"] or 1) <= 1.01 and (f4["costo"] or 1) <= 1.01:
            return {"variante": V, "p_max": 24, "costo": 1,
                    "motivo": "p hilos en p nucleos fisicos, Slurm cobra p CPUs, cabe hasta p = 24"}
    for V in ("C", "B"):
        f4, f24 = ver(V, 4), ver(V, 24)
        if pasa(f4) and (f4["costo"] or 1) <= 1.01 and not pasa(f24):
            por_que = "no se probo" if f24 is None else "%s (%s)" % (f24["veredicto"], f24["motivo"][:60])
            return {"variante": V, "p_max": 12, "costo": 1,
                    "motivo": "funciona con p = 4, pero con p = 24 %s; se mide hasta p = 12" % por_que}
    f4 = ver("D", 4)
    if pasa(f4):
        return {"variante": "D", "p_max": LIMITE_QOS // 2, "costo": 2,
                "motivo": "solo funciona pidiendo 2 hilos de hardware por nucleo: cuesta el doble y cabe hasta p = 12"}
    return None


def imprimir(filas, rec):
    print("Diagnostico de OpenMP (Tarea A.2): p hilos en p nucleos fisicos?")
    print()
    for V in "ABCD":
        print("  %s : %s" % (V, DESCRIPCION[V]))
    print()
    cab = "%-4s %-8s %-6s %-9s %-9s %-7s %-6s %-9s %-10s %s" % (
        "Var", "Job", "CPUs", "Mascara", "Hilos", "VASP", "Iter", "t/iter", "Veredicto", "Detalle")
    print(cab)
    print("-" * len(cab))
    for f in filas:
        print("%-4s %-8s %-6s %-9s %-9s %-7s %-6s %-9s %-10s %s" % (
            "%s%d" % (f["variante"], f["p"]), f["jobid"] if f["jobid"] != "RECHAZADO" else "-",
            f.get("alloc") or "-",
            "-" if f.get("nuc_mascara") is None else "%d/%d" % (f["nuc_mascara"], f["p"]),
            "-" if not f.get("muestras") else "%d/%d" % (f["nuc_hilos"], f["p"]),
            f.get("hilos_vasp") or "-", f.get("n_iter") or "-",
            "-" if f.get("t_iter_s") is None else "%.1f s" % f["t_iter_s"],
            f["veredicto"], f["motivo"]))
    print()
    print("CPUs = CPUs que Slurm asigno (y cobra). Mascara = nucleos fisicos permitidos al proceso.")
    print("Hilos = nucleos fisicos donde corrieron de verdad los hilos de VASP (ps -L durante la corrida).")
    print()
    pend = [f for f in filas if f["veredicto"] in ("EN COLA", "CORRIENDO")]
    if pend:
        print("Todavia hay corridas en cola o corriendo: %s. Espera y repite."
              % ", ".join("%s%d" % (f["variante"], f["p"]) for f in pend))
    elif rec:
        print("RECOMENDACION: variante %s  (%s)" % (rec["variante"], DESCRIPCION[rec["variante"]]))
        print("               %s." % rec["motivo"])
        print("Escrita en datos/parteA2_variante.txt. Para guardar la evidencia en Git:")
        print("  python3 analisis/parteA2_diagnostico.py --guardar")
    else:
        print("Ninguna variante funciono: no envies el barrido y mandame una captura.")


def guardar(filas, rec):
    dd = os.path.join(REPO, "datos")
    os.path.isdir(dd) or os.makedirs(dd)
    cols = ["variante", "p", "jobid", "alloc", "costo", "nuc_mascara", "nuc_hilos", "muestras", "hilos_vasp",
            "ranks", "n_iter", "t_iter_s", "socket", "veredicto", "motivo"]
    with io.open(os.path.join(dd, "parteA2_diagnostico.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for f in filas:
            w.writerow(["" if f.get(c) is None else (("%.3f" % f[c]) if isinstance(f[c], float) else f[c]) for c in cols])
    md = [u"# Tarea A.2 - diagnostico de OpenMP en Granado", u"",
          u"Corridas cortas del caso base (NELM = 4). Pregunta: los p hilos de VASP corren en p nucleos fisicos?",
          u"", u"| Variante | Peticion a Slurm | p | CPUs cobradas | Nucleos en la mascara | Nucleos usados por los hilos | Veredicto |",
          u"|---|---|---|---|---|---|---|"]
    for f in filas:
        md.append(u"| %s | `%s` | %d | %s | %s | %s | %s: %s |" % (
            f["variante"], DESCRIPCION.get(f["variante"], "?"), f["p"], f.get("alloc") or "-",
            "-" if f.get("nuc_mascara") is None else f["nuc_mascara"],
            "-" if not f.get("muestras") else f["nuc_hilos"], f["veredicto"], f["motivo"]))
    if rec:
        md += [u"", u"**Variante elegida: %s** (`%s`): %s." % (rec["variante"], DESCRIPCION[rec["variante"]], rec["motivo"])]
    with io.open(os.path.join(dd, "parteA2_diagnostico.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md) + u"\n")
    if rec:
        with io.open(SALIDA, "w", encoding="utf-8") as fh:
            fh.write(u"variante=%s\np_max=%d\ncosto=%d\nmotivo=%s\n" % (rec["variante"], rec["p_max"], rec["costo"], rec["motivo"]))
    elif os.path.exists(SALIDA):
        os.remove(SALIDA)


def guardar_evidencia(filas):
    base = os.path.join(REPO, "datos", "parteA2_diag")
    for f in filas:
        d = f.get("dir")
        if not d or not os.path.isdir(d):
            continue
        dst = os.path.join(base, os.path.basename(d))
        os.path.isdir(dst) or os.makedirs(dst)
        for a in ("tiempo.txt", "entorno.txt", "afinidad.txt", "hilos.txt", "hilos_crudo.txt", "slurm_job.txt",
                  "carga.txt", "variables.txt"):
            if os.path.isfile(os.path.join(d, a)):
                shutil.copy(os.path.join(d, a), os.path.join(dst, a))
        err = leer(os.path.join(d, "vasp.err")).splitlines()
        if err:
            with io.open(os.path.join(dst, "vasp_err_inicio.txt"), "w", encoding="utf-8") as fh:
                fh.write(u"\n".join(err[:60]) + u"\n")
    out = git("add", "datos")
    if out.startswith("ERROR"):
        print("AVISO: git add fallo: " + out)
        return
    if git("diff", "--cached", "--quiet").startswith("ERROR"):
        out = git("commit", "-q", "-m", "Parte A.2: diagnostico de OpenMP y variante elegida")
        print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
    else:
        print("Nada nuevo que guardar en Git.")


def main():
    ap = argparse.ArgumentParser(description="Resultado del diagnostico de OpenMP (A.2)")
    ap.add_argument("--guardar", action="store_true", help="guarda la evidencia en datos/ y hace commit")
    a = ap.parse_args()
    if not os.path.isfile(REG):
        print("Todavia no se envio el diagnostico: primero  bash slurm/parteA2_diagnostico.sh")
        return
    filas = evaluar()
    pend = [f for f in filas if f["veredicto"] in ("EN COLA", "CORRIENDO")]
    rec = None if pend else recomendar(filas)
    imprimir(filas, rec)
    guardar(filas, rec)
    if a.guardar:
        if pend:
            print("No se guarda nada mientras haya corridas pendientes.")
        else:
            print()
            guardar_evidencia(filas)


if __name__ == "__main__":
    main()
