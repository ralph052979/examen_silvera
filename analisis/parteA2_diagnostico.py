#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteA2_diagnostico.py - Tarea A.2, ronda 1 del diagnostico de OpenMP: como pedir
                         p nucleos FISICOS para p hilos (variantes A, B, C y D).

Lee las 6 corridas cortas de logs/jobs_A2d.txt (NELM = 4). Para cada una informa la
mascara que dio srun (nucleos fisicos), cuantas CPUs cobro Slurm, si VASP termino las
4 iteraciones y cuanto tardo cada iteracion frente a las corridas con p = 1 de la A.1.
Las corridas que Slurm corto por tiempo (TIMEOUT) no dejan tiempo.txt: su mascara se
lee directamente de vasp.err.

Uso (desde la carpeta del repositorio):
  python3 analisis/parteA2_diagnostico.py             # tabla y conclusiones
  python3 analisis/parteA2_diagnostico.py --guardar   # ademas guarda la evidencia en datos/
                                                      # (y datos/parteA2_geometria.txt) y hace commit

Solo usa la biblioteca estandar de Python. Todo lo que imprime es ASCII.
"""
from __future__ import print_function
import argparse
import csv
import glob
import io
import os
import re
import shutil
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REG = os.path.join(REPO, "logs", "jobs_A2d.txt")
GEOMETRIA = os.path.join(REPO, "datos", "parteA2_geometria.txt")
VIEJO = os.path.join(REPO, "datos", "parteA2_variante.txt")
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


# ----------------------------------------------------------------------------- topologia
def topologia():
    """cpu -> (socket, nucleo). Se usa topologia.txt de cualquier corrida que lo tenga
    (las de la ronda 2 en adelante); si no hay, la de los nodos de Granado:
    2 sockets x 18 nucleos x 2 hilos; el CPU c y el c + 36 son el mismo nucleo fisico."""
    for f in sorted(glob.glob(os.path.join(REPO, "corridas", "*", "topologia.txt"))):
        t = {}
        for l in leer(f).splitlines():
            p = l.split()
            if len(p) == 3 and p[0].isdigit():
                t[int(p[0])] = (p[1], p[2])
        if len(t) >= 2:
            return t, "topologia.txt de %s" % os.path.basename(os.path.dirname(f))
    return dict((c, (str((c % 36) // 18), str(c % 18))) for c in range(72)), "regla de Granado (c y c+36 = mismo nucleo)"


def mascaras(err, topo):
    """lineas de srun --cpu-bind=verbose -> lista de (tarea, cpus logicos, nucleos fisicos)"""
    out = []
    vistos = set()
    for m in re.finditer(r"task\s+(\d+)\s+\d+\s+\[\d+\]:\s+mask\s+0x([0-9a-fA-F]+)", err):
        tk = int(m.group(1))
        if tk in vistos:
            continue
        vistos.add(tk)
        v = int(m.group(2), 16)
        cpus = [b for b in range(v.bit_length()) if (v >> b) & 1]
        nuc = set(topo.get(c, ("?", str(c))) for c in cpus)
        out.append((tk, cpus, nuc))
    return out


def leer_outcar(ruta):
    t = leer(ruta)
    o = {}
    m = re.search(r"running\s+(\d+)\s+mpi-ranks,\s+with\s+(\d+)\s+threads/rank", t)
    o["ranks"] = int(m.group(1)) if m else None
    o["hilos_vasp"] = int(m.group(2)) if m else None
    pares = [(float(a), float(b)) for a, b in re.findall(r"LOOP:\s+cpu time\s+([\d.]+)\s*:\s*real time\s+([\d.]+)", t)]
    o["n_iter"] = len(pares)
    o["t_iter_s"] = sum(b for a, b in pares) / len(pares) if pares else None
    o["cpu_real"] = (sum(a for a, b in pares) / sum(b for a, b in pares)) if pares and sum(b for a, b in pares) > 0 else None
    o["loops"] = pares
    return o


def referencia_a1():
    """t por iteracion de las primeras 4 iteraciones de las corridas con p = 1 de la A.1"""
    vals = []
    for d in sorted(glob.glob(os.path.join(REPO, "corridas", "ex_A1_p1_r*_*"))):
        if clave_valor(os.path.join(d, "tiempo.txt")).get("estado") != "OK":
            continue
        o = leer_outcar(os.path.join(d, "OUTCAR"))
        vals += [b for a, b in o["loops"][:4]]
    return (sum(vals) / len(vals), len(vals)) if vals else (None, 0)


# ----------------------------------------------------------------------------- evaluacion
def evaluar():
    topo, fuente = topologia()
    filas = []
    for l in leer(REG).splitlines():
        p = l.split()
        if len(p) < 3:
            continue
        V, P, jid = p[0], int(p[1]), p[2]
        f = {"variante": V, "p": P, "jobid": jid, "veredicto": "", "motivo": "", "alloc": None,
             "nuc_mascara": None, "cpus_mascara": None, "n_iter": None, "t_iter_s": None, "cpu_real": None,
             "estado": None, "slurm": None}
        if jid == "RECHAZADO":
            f["veredicto"] = "RECHAZADO"
            f["motivo"] = " ".join(p[3:])[:110]
            filas.append(f)
            continue
        d = os.path.join(REPO, "corridas", "ex_A2d_%s%d_%s" % (V, P, jid))
        f["dir"] = d
        st, mot = slurm(jid)
        f["slurm"] = st or "?"
        if st in ("PENDING", "CONFIGURING", "RUNNING", "COMPLETING") or (st is None and not os.path.isdir(d)):
            f["veredicto"] = "EN CURSO"
            f["motivo"] = mot or ""
            filas.append(f)
            continue
        tie = clave_valor(os.path.join(d, "tiempo.txt"))
        f["estado"] = tie.get("estado")
        m = re.search(r"NumCPUs=(\d+)", leer(os.path.join(d, "slurm_job.txt")))
        f["alloc"] = int(m.group(1)) if m else None
        ms = mascaras(leer(os.path.join(d, "vasp.err")), topo)
        if ms:
            f["nuc_mascara"] = len(ms[0][2])
            f["cpus_mascara"] = len(ms[0][1])
            f["cpus_lista"] = ms[0][1]
        f.update(dict((k, v) for k, v in leer_outcar(os.path.join(d, "OUTCAR")).items() if k != "loops"))
        problemas = []
        if f["nuc_mascara"] is not None and f["nuc_mascara"] < P:
            problemas.append("mascara con %d nucleos fisicos de %d" % (f["nuc_mascara"], P))
        if f["estado"] is None:
            problemas.append("Slurm la corto (%s): no termino las 4 iteraciones en 20 min" % (st or "sin estado"))
        elif f["estado"] != "OK":
            problemas.append(f["estado"])
        if problemas:
            f["veredicto"] = "NO PASA"
            f["motivo"] = "; ".join(problemas)
        else:
            f["veredicto"] = "TERMINO"
            f["motivo"] = "mascara de %d nucleos fisicos, 4 iteraciones completas" % P
        filas.append(f)
    return filas, fuente


def imprimir(filas, fuente, ref):
    print("Diagnostico de OpenMP, ronda 1 (Tarea A.2): como pedir p nucleos fisicos para p hilos")
    print()
    for V in "ABCD":
        print("  %s : %s" % (V, DESCRIPCION[V]))
    print()
    cab = "%-4s %-8s %-10s %-5s %-14s %-5s %-9s %-8s %-9s %s" % (
        "Var", "Job", "Slurm", "CPUs", "Mascara", "Iter", "t/iter", "CPU/real", "Veredicto", "Detalle")
    print(cab)
    print("-" * len(cab))
    for f in filas:
        masc = "-" if f.get("nuc_mascara") is None else "%d nuc/%d log" % (f["nuc_mascara"], f["cpus_mascara"])
        print("%-4s %-8s %-10s %-5s %-14s %-5s %-9s %-8s %-9s %s" % (
            "%s%d" % (f["variante"], f["p"]), f["jobid"] if f["jobid"] != "RECHAZADO" else "-",
            (f.get("slurm") or "-")[:10], f.get("alloc") or "-", masc, f.get("n_iter") or "-",
            "-" if f.get("t_iter_s") is None else "%.1f s" % f["t_iter_s"],
            "-" if f.get("cpu_real") is None else "%.2f" % f["cpu_real"],
            f["veredicto"], f["motivo"]))
    print()
    print("Mascara = nucleos fisicos / CPUs logicos que srun permitio al proceso (topologia: %s)." % fuente)
    print("CPU/real = tiempo de CPU / tiempo real de las iteraciones (con p hilos bien usados vale ~p).")
    if ref[0]:
        print("Referencia: con 1 proceso y 1 hilo (A.1, p = 1) cada una de las primeras 4 iteraciones"
              " tarda %.1f s (promedio de %d)." % ref)
    print("Las corridas que Slurm corto no dejan OSZICAR ni OUTCAR completos: VASP compilado con gfortran")
    print("los escribe en bloques y solo los vacia al terminar.")
    print()
    pend = [f for f in filas if f["veredicto"] == "EN CURSO"]
    if pend:
        print("Todavia hay corridas en cola o corriendo: %s. Espera y repite."
              % ", ".join("%s%d" % (f["variante"], f["p"]) for f in pend))
        return
    print("CONCLUSIONES DE LA RONDA 1:")
    for f in filas:
        if f["veredicto"] == "TERMINO" and f.get("t_iter_s") and ref[0]:
            print("  * %s%d termino, pero cada iteracion tardo %.1f s: %.2f veces el tiempo con 1 hilo (%.1f s),"
                  % (f["variante"], f["p"], f["t_iter_s"], f["t_iter_s"] / ref[0], ref[0]))
            if f.get("cpu_real"):
                print("    con un uso de CPU de %.2f para %d hilos: hay mas hilos trabajando que los de OpenMP."
                      % (f["cpu_real"], f["p"]))
    cortadas = ["%s%d" % (f["variante"], f["p"]) for f in filas if f["veredicto"] == "NO PASA" and f.get("estado") is None]
    if cortadas:
        print("  * %s no terminaron 4 iteraciones en 20 min (con 1 hilo bastan unos 7 min)." % ", ".join(cortadas))
    mal = ["%s%d (%d de %d)" % (f["variante"], f["p"], f["nuc_mascara"], f["p"]) for f in filas
           if f.get("nuc_mascara") is not None and f["nuc_mascara"] < f["p"]]
    if mal:
        print("  * Mascara con menos nucleos fisicos que hilos: %s." % ", ".join(mal))
    print("  * Ninguna variante sirve tal como esta: hace falta la ronda 2 (bash slurm/parteA2_hilos.sh),")
    print("    que prueba la causa (los hilos propios de OpenBLAS y la espera activa de OpenMP).")
    g = geometria(filas)
    print("  * Para la ronda 2 se usa la variante B%s." % (
        "" if g["paso_tpc"] == 0 else ", y ademas srun --threads-per-core=1 (la mascara de B no dio p nucleos)"))


def geometria(filas):
    def ver(V, P):
        for f in filas:
            if f["variante"] == V and f["p"] == P:
                return f
        return None
    b4, b24 = ver("B", 4), ver("B", 24)
    n4 = b4.get("nuc_mascara") if b4 else None
    n24 = b24.get("nuc_mascara") if b24 else None
    tpc = 0 if (n4 == 4 and n24 == 24) else 1
    return {"b4_nucleos": n4, "b24_nucleos": n24, "paso_tpc": tpc}


def guardar(filas, fuente, ref):
    dd = os.path.join(REPO, "datos")
    os.path.isdir(dd) or os.makedirs(dd)
    cols = ["variante", "p", "jobid", "slurm", "alloc", "nuc_mascara", "cpus_mascara", "hilos_vasp", "n_iter",
            "t_iter_s", "cpu_real", "estado", "veredicto", "motivo"]
    with io.open(os.path.join(dd, "parteA2_diagnostico.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for f in filas:
            w.writerow(["" if f.get(c) is None else (("%.3f" % f[c]) if isinstance(f[c], float) else f[c]) for c in cols])
    md = [u"# Tarea A.2 - diagnostico de OpenMP, ronda 1 (como pedir los nucleos)", u"",
          u"Corridas cortas del caso base (NELM = 4), enviadas a la vez. Topologia: %s." % fuente, u"",
          u"| Variante | Peticion a Slurm | p | CPUs cobradas | Mascara (nucleos fisicos / CPUs logicos) | Iteraciones | t/iter (s) | CPU/real | Veredicto |",
          u"|---|---|---|---|---|---|---|---|---|"]
    for f in filas:
        md.append(u"| %s | `%s` | %d | %s | %s | %s | %s | %s | %s: %s |" % (
            f["variante"], DESCRIPCION.get(f["variante"], "?"), f["p"], f.get("alloc") or "-",
            "-" if f.get("nuc_mascara") is None else "%d / %d" % (f["nuc_mascara"], f["cpus_mascara"]),
            f.get("n_iter") or "-", "-" if f.get("t_iter_s") is None else "%.1f" % f["t_iter_s"],
            "-" if f.get("cpu_real") is None else "%.2f" % f["cpu_real"], f["veredicto"], f["motivo"]))
    if ref[0]:
        md += [u"", u"Referencia: 1 proceso con 1 hilo (A.1, p = 1), primeras 4 iteraciones: %.1f s por iteracion." % ref[0]]
    md += [u"", u"Las corridas que Slurm corto por tiempo no dejan OSZICAR ni OUTCAR: VASP compilado con gfortran "
                u"los vacia al terminar. Su mascara se leyo de `vasp.err` (lineas de `srun --cpu-bind=verbose`)."]
    with io.open(os.path.join(dd, "parteA2_diagnostico.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md) + u"\n")
    g = geometria(filas)
    with io.open(GEOMETRIA, "w", encoding="utf-8") as fh:
        fh.write(u"b4_nucleos=%s\nb24_nucleos=%s\npaso_tpc=%d\n" % (g["b4_nucleos"], g["b24_nucleos"], g["paso_tpc"]))


def guardar_evidencia(filas):
    base = os.path.join(REPO, "datos", "parteA2_diag")
    for f in filas:
        d = f.get("dir")
        if not d or not os.path.isdir(d):
            continue
        dst = os.path.join(base, os.path.basename(d))
        os.path.isdir(dst) or os.makedirs(dst)
        for a in ("tiempo.txt", "entorno.txt", "afinidad.txt", "hilos.txt", "hilos_crudo.txt", "slurm_job.txt",
                  "carga.txt", "variables.txt", "OSZICAR"):
            if os.path.isfile(os.path.join(d, a)):
                shutil.copy(os.path.join(d, a), os.path.join(dst, a))
        err = leer(os.path.join(d, "vasp.err")).splitlines()
        if err:
            with io.open(os.path.join(dst, "vasp_err_inicio.txt"), "w", encoding="utf-8") as fh:
                fh.write(u"\n".join(err[:60]) + u"\n")
        loops = [l for l in leer(os.path.join(d, "OUTCAR")).splitlines() if "LOOP:" in l or "running" in l]
        if loops:
            with io.open(os.path.join(dst, "outcar_loop.txt"), "w", encoding="utf-8") as fh:
                fh.write(u"\n".join(loops) + u"\n")
    if os.path.exists(VIEJO):
        git("rm", "-q", "--cached", "--ignore-unmatch", "datos/parteA2_variante.txt")
        os.remove(VIEJO)
    out = git("add", "datos")
    if out.startswith("ERROR"):
        print("AVISO: git add fallo: " + out)
        return
    if git("diff", "--cached", "--quiet").startswith("ERROR"):
        out = git("commit", "-q", "-m", "Parte A.2: diagnostico de OpenMP, ronda 1 (mascaras y tiempos)")
        print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
    else:
        print("Nada nuevo que guardar en Git.")


def main():
    ap = argparse.ArgumentParser(description="Diagnostico de OpenMP, ronda 1 (A.2)")
    ap.add_argument("--guardar", action="store_true", help="guarda la evidencia en datos/ y hace commit")
    a = ap.parse_args()
    if not os.path.isfile(REG):
        print("No hay registro de la ronda 1 (logs/jobs_A2d.txt).")
        return
    filas, fuente = evaluar()
    ref = referencia_a1()
    imprimir(filas, fuente, ref)
    pend = [f for f in filas if f["veredicto"] == "EN CURSO"]
    if a.guardar:
        print()
        if pend:
            print("No se guarda nada mientras haya corridas en cola o corriendo.")
        else:
            guardar(filas, fuente, ref)
            guardar_evidencia(filas)


if __name__ == "__main__":
    main()
