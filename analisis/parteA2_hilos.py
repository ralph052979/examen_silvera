#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteA2_hilos.py - Tarea A.2, ronda 2 del diagnostico de OpenMP: hilos de OpenBLAS y
                   politica de espera de OpenMP (diseno 2 x 2 con p = 4, y p = 24).

Para cada corrida de logs/jobs_A2h.txt informa: hilos ACTIVOS de VASP (>= 50 % de un
CPU) frente a los nucleos que tiene, uso total de CPU, tiempo por iteracion SCF y la
aceleracion por iteracion frente a 1 hilo (primeras 4 iteraciones de la A.1, p = 1).
Elige el entorno para el barrido (E1 si funciona; si no, E3), estima su costo y lo
escribe en datos/parteA2_config.txt.

Uso (desde la carpeta del repositorio):
  python3 analisis/parteA2_hilos.py             # tabla, explicacion y eleccion
  python3 analisis/parteA2_hilos.py --guardar   # ademas guarda la evidencia en datos/ y hace commit

Solo usa la biblioteca estandar de Python. Todo lo que imprime es ASCII.
"""
from __future__ import print_function
import argparse
import csv
import glob
import io
import math
import os
import re
import shutil
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REG = os.path.join(REPO, "logs", "jobs_A2h.txt")
GEO = os.path.join(REPO, "datos", "parteA2_geometria.txt")
CONFIG = os.path.join(REPO, "datos", "parteA2_config.txt")
P_BARRIDO = [2, 4, 8, 12, 16, 24]
REPETICIONES = 3
ENTORNO = {
    "E0": ("de la ronda 1 (control)", "", ""),
    "E1": ("OPENBLAS_NUM_THREADS=1", "1", ""),
    "E2": ("OMP_WAIT_POLICY=PASSIVE", "", "PASSIVE"),
    "E3": ("OPENBLAS_NUM_THREADS=1 + OMP_WAIT_POLICY=PASSIVE", "1", "PASSIVE"),
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


def num(x, tipo=float):
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


def estado_slurm(jobid):
    for cmd in (["squeue", "-h", "-j", str(jobid), "-o", "%T"], ["sacct", "-n", "-X", "-j", str(jobid), "-o", "State%20"]):
        try:
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT).decode().split()
        except (OSError, subprocess.CalledProcessError):
            continue
        if out:
            return out[0]
    return None


def loops(ruta):
    t = leer(ruta)
    return [(float(a), float(b)) for a, b in re.findall(r"LOOP:\s+cpu time\s+([\d.]+)\s*:\s*real time\s+([\d.]+)", t)]


def hilos_vasp(ruta):
    m = re.search(r"running\s+(\d+)\s+mpi-ranks,\s+with\s+(\d+)\s+threads/rank", leer(ruta))
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)


def referencia_a1():
    """(t por iteracion de las primeras 4 iteraciones, t de pared medio) de la A.1 con p = 1"""
    its, walls = [], []
    for d in sorted(glob.glob(os.path.join(REPO, "corridas", "ex_A1_p1_r*_*"))):
        tie = clave_valor(os.path.join(d, "tiempo.txt"))
        if tie.get("estado") != "OK":
            continue
        its += [b for a, b in loops(os.path.join(d, "OUTCAR"))[:4]]
        w = num(tie.get("wall_s"))
        if w:
            walls.append(w)
    t_it = sum(its) / len(its) if its else None
    t1 = sum(walls) / len(walls) if walls else None
    return t_it, t1, len(walls)


def resumen_pila(ruta):
    t = leer(ruta)
    if not t:
        return None
    if "no esta instalado" in t:
        return "gdb no esta instalado en el nodo"
    bloques = re.split(r"\nThread \d+ ", "\n" + t)[1:]
    if not bloques:
        return "sin pilas (%s)" % t.strip().splitlines()[-1][:70] if t.strip() else "vacio"
    n_ob = sum(1 for b in bloques if "openblas" in b.lower() or "blas_thread_server" in b)
    n_omp = sum(1 for b in bloques if "gomp" in b.lower() and "openblas" not in b.lower())
    return "%d hilos: %d dentro de OpenBLAS, %d de OpenMP (libgomp), %d otros" % (
        len(bloques), n_ob, n_omp, len(bloques) - n_ob - n_omp)


def evaluar(t_ref):
    filas = []
    for l in leer(REG).splitlines():
        p = l.split()
        if len(p) < 3:
            continue
        E, P, jid = p[0], int(p[1]), p[2]
        f = {"E": E, "p": P, "jobid": jid, "veredicto": "", "motivo": "", "t_iter_s": None, "S_iter": None,
             "n_iter": None, "activos": None, "uso": None, "nuc_mascara": None, "estado": None}
        if jid == "RECHAZADO":
            f["veredicto"] = "RECHAZADO"
            f["motivo"] = " ".join(p[3:])[:100]
            filas.append(f)
            continue
        d = os.path.join(REPO, "corridas", "ex_A2h_%s_p%d_%s" % (E, P, jid))
        f["dir"] = d
        tie = clave_valor(os.path.join(d, "tiempo.txt"))
        if not tie.get("estado"):
            st = estado_slurm(jid)
            if st in (None, "PENDING", "RUNNING", "CONFIGURING", "COMPLETING"):
                f["veredicto"] = "EN CURSO"
                f["motivo"] = st or "en cola"
            else:
                f["veredicto"] = "NO PASA"
                f["motivo"] = "el job termino como %s sin resultados" % st
            filas.append(f)
            continue
        f["estado"] = tie["estado"]
        if f["estado"].startswith("OMITIDA"):
            f["veredicto"] = "OMITIDA"
            f["motivo"] = "se omitio: " + f["estado"][8:].replace("_", " ")[:80]
            filas.append(f)
            continue
        ent = clave_valor(os.path.join(d, "entorno.txt"))
        f["openblas_omp"] = ent.get("openblas_con_openmp")
        f["openblas_ver"] = ent.get("openblas_version")
        f["fftlib"] = ent.get("cadena_fftlib_en_vasp")
        f["nodo"] = (ent.get("nodo") or "-").split(".")[0]
        f["nuc_mascara"] = num(tie.get("nucleos_fisicos_distintos"), int) or None
        f["intervalos"] = num(tie.get("hilos_intervalos"), int) or 0
        f["activos"] = num(tie.get("hilos_activos_mediana")) if f["intervalos"] else None
        f["uso"] = num(tie.get("uso_cpu_mediana")) if f["intervalos"] else None
        f["hilos_totales"] = num(tie.get("hilos_totales"), int)
        f["alloc"] = num(tie.get("alloc_cpus"), int)
        f["carga_ajena"] = tie.get("carga_ajena_max", "-")
        lp = loops(os.path.join(d, "OUTCAR"))
        ranks, hv = hilos_vasp(os.path.join(d, "OUTCAR"))
        f["hilos_vasp"] = hv
        if lp:
            f["n_iter"] = len(lp)
            f["t_iter_s"] = sum(b for a, b in lp) / len(lp)
            f["cpu_real"] = sum(a for a, b in lp) / sum(b for a, b in lp)
        else:
            prog = [num(x.split()[2]) for x in leer(os.path.join(d, "progreso.txt")).splitlines() if len(x.split()) == 3]
            prog = [x for x in prog if x is not None]
            f["n_iter"] = len(prog)
            if len(prog) >= 2:
                f["t_iter_s"] = (prog[-1] - prog[0]) / (len(prog) - 1)
            elif tie.get("detenida_por") == "lento":
                f["t_iter_s"] = None
                f["sin_iter_en"] = num(re.sub(r"\D", "", f["estado"].split("en_")[-1]) or None)
        if f["t_iter_s"] and t_ref:
            f["S_iter"] = t_ref / f["t_iter_s"]
        f["pila"] = resumen_pila(os.path.join(d, "pila_hilos.txt"))
        if f["estado"] == "OK" and f.get("S_iter") is not None and f["S_iter"] < 1.0:
            f["veredicto"] = "LENTA"
            f["motivo"] = "mas lenta que con 1 hilo (%g hilos activos para %d nucleos)" % (f["activos"] or 0, P)
        elif f["estado"] == "OK":
            f["veredicto"] = "PASA"
            f["motivo"] = "%g hilos activos para %d nucleos" % (f["activos"], P) if f["activos"] is not None else "OK"
        else:
            f["veredicto"] = "NO PASA"
            f["motivo"] = f["estado"].replace("FALLA_", "").replace("_", " ")[:90]
            if f["estado"].startswith("FALLA_srun"):
                err = [x.strip() for x in leer(os.path.join(d, "vasp.err")).splitlines() if "error" in x.lower()]
                if err:
                    f["motivo"] += ": " + err[0][:90]
        filas.append(f)
    return filas


def ver(filas, E, P):
    for f in filas:
        if f["E"] == E and f["p"] == P:
            return f
    return None


def elegir(filas):
    for E in ("E1", "E3"):
        a, b = ver(filas, E, 4), ver(filas, E, 24)
        if a and b and a["veredicto"] == "PASA" and b["veredicto"] == "PASA":
            return E
    return None


def proyectar(filas, E, t_ref, t1):
    """Amdahl por minimos cuadrados con S por iteracion en p = 4 y 24; t(p) = t1 / S(p)"""
    pts = [(f["p"], f["S_iter"]) for f in (ver(filas, E, 4), ver(filas, E, 24)) if f and f.get("S_iter")]
    if not pts or not t1:
        return None
    sxx = sum((1 - 1.0 / p) ** 2 for p, s in pts)
    sxy = sum((1 - 1.0 / p) * (1.0 / s - 1.0 / p) for p, s in pts)
    fa = min(1.0, max(0.0, sxy / sxx))
    res = {"f": fa, "t": {}, "lim": {}}
    costo = 0.0
    for p in P_BARRIDO:
        S = 1.0 / (fa + (1 - fa) / p)
        t = t1 / S
        res["t"][p] = t
        costo += REPETICIONES * p * t / 3600.0
        lim = max(20 * 60, 1.5 * t + 600)
        lim = int(math.ceil(lim / 300.0) * 300)
        lim = min(lim, 4 * 3600 + 55 * 60)
        res["lim"][p] = "%02d:%02d:00" % (lim // 3600, (lim % 3600) // 60)
    res["costo"] = costo
    return res


def imprimir(filas, t_ref, t1, n1, eleccion, proy, geo):
    print("Diagnostico de OpenMP, ronda 2 (Tarea A.2): hilos de OpenBLAS y espera de OpenMP")
    print()
    for E in ("E0", "E1", "E2", "E3"):
        print("  %s : %s" % (E, ENTORNO[E][0]))
    print("  (1 proceso, variante B: --threads-per-core=1%s)" % (" + srun --threads-per-core=1" if geo.get("paso_tpc") == "1" else ""))
    print()
    cab = "%-4s %-3s %-8s %-8s %-9s %-6s %-5s %-9s %-7s %-9s %s" % (
        "Ent", "p", "Job", "Mascara", "Activos", "UsoCPU", "Iter", "t/iter", "S_iter", "Veredicto", "Detalle")
    print(cab)
    print("-" * len(cab))
    for f in filas:
        if f.get("t_iter_s"):
            ti = "%.1f s" % f["t_iter_s"]
        elif f.get("sin_iter_en"):
            ti = "> %d s" % f["sin_iter_en"]
        else:
            ti = "-"
        print("%-4s %-3s %-8s %-8s %-9s %-6s %-5s %-9s %-7s %-9s %s" % (
            f["E"], f["p"], f["jobid"] if f["jobid"] != "RECHAZADO" else "-",
            "-" if f.get("nuc_mascara") is None else "%d/%d" % (f["nuc_mascara"], f["p"]),
            "-" if f.get("activos") is None else "%g/%d" % (f["activos"], f["p"]),
            "-" if f.get("uso") is None else "%.2f" % f["uso"],
            "-" if f.get("n_iter") is None else f["n_iter"], ti,
            "-" if f.get("S_iter") is None else "%.2f" % f["S_iter"], f["veredicto"], f["motivo"]))
    print()
    print("Mascara = nucleos fisicos del proceso / hilos. Activos = mediana de hilos con >= 50 % de un CPU.")
    print("UsoCPU = CPUs que usa el proceso (mediana). S_iter = t por iteracion con 1 hilo / t por iteracion.")
    print("PASA = mascara de p nucleos, sin hilos de mas y no mas lenta que con 1 hilo.")
    if t_ref:
        print("Referencia con 1 hilo (A.1, p = 1, primeras 4 iteraciones): %.1f s por iteracion; corrida completa: %.0f s (%d corridas)."
              % (t_ref, t1 or 0, n1))
    ob = [f for f in filas if f.get("openblas_omp")]
    if ob:
        g = ob[0]
        print("OpenBLAS: %s; %s." % (g.get("openblas_ver") or "version ?",
              "usa los hilos de OpenMP (enlaza libgomp)" if g["openblas_omp"] == "si" else
              "usa hilos propios (pthreads, no enlaza libgomp)" if g["openblas_omp"] == "no" else "no se pudo saber si usa OpenMP"))
    pilas = [f for f in filas if f.get("pila")]
    for f in pilas:
        print("Pila de los hilos (gdb) en %s p=%d: %s." % (f["E"], f["p"], f["pila"]))
    print()
    pend = [f for f in filas if f["veredicto"] == "EN CURSO"]
    if pend:
        print("Todavia hay corridas en cola o corriendo: %s. Espera y repite."
              % ", ".join("%s p=%d" % (f["E"], f["p"]) for f in pend))
        return
    e0, e1 = ver(filas, "E0", 4), ver(filas, "E1", 4)
    if e0 and e1 and e0.get("activos") is not None and e1.get("activos") is not None:
        print("Con el entorno de la ronda 1 (E0) trabajan %g hilos para 4 nucleos; con OPENBLAS_NUM_THREADS=1 (E1), %g."
              % (e0["activos"], e1["activos"]))
    c1, c3 = ver(filas, "E1", 24), ver(filas, "E3", 24)
    if c1 and c3 and c1.get("t_iter_s") and c3.get("t_iter_s"):
        print("Con p = 24, la espera pasiva (E3) cambia el tiempo por iteracion en %+.0f %% respecto de E1."
              % (100.0 * (c3["t_iter_s"] - c1["t_iter_s"]) / c1["t_iter_s"]))
    if eleccion:
        a, b = ver(filas, eleccion, 4), ver(filas, eleccion, 24)
        print("ELECCION: %s (%s): pasa con p = 4 y p = 24." % (eleccion, ENTORNO[eleccion][0]))
        if a.get("S_iter") and b.get("S_iter"):
            print("  Aceleracion por iteracion: S(4) = %.2f, S(24) = %.2f." % (a["S_iter"], b["S_iter"]))
        if proy:
            print("  Amdahl con esos dos puntos: f = %.3f. Barrido proyectado (3 repeticiones de p = %s):"
                  % (proy["f"], ", ".join(str(p) for p in P_BARRIDO)))
            print("  " + "  ".join("p=%d: %.0f s" % (p, proy["t"][p]) for p in P_BARRIDO))
            print("  costo estimado %.0f h-nucleo." % proy["costo"])
            if proy["costo"] > 40:
                print("  AVISO: es caro. Mandame una captura antes de enviar el barrido.")
        print("Escrito en datos/parteA2_config.txt. Para guardar la evidencia en Git:")
        print("  python3 analisis/parteA2_hilos.py --guardar")
    else:
        print("Ningun entorno paso con p = 4 y p = 24: no envies el barrido y mandame una captura.")


def guardar(filas, t_ref, t1, eleccion, proy, geo):
    dd = os.path.join(REPO, "datos")
    os.path.isdir(dd) or os.makedirs(dd)
    cols = ["E", "p", "jobid", "nodo", "alloc", "nuc_mascara", "hilos_vasp", "hilos_totales", "activos", "uso",
            "n_iter", "t_iter_s", "cpu_real", "S_iter", "carga_ajena", "openblas_omp", "estado", "veredicto", "motivo"]
    with io.open(os.path.join(dd, "parteA2_hilos.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for f in filas:
            w.writerow(["" if f.get(c) is None else (("%.3f" % f[c]) if isinstance(f[c], float) else f[c]) for c in cols])
    md = [u"# Tarea A.2 - diagnostico de OpenMP, ronda 2 (hilos de OpenBLAS y espera de OpenMP)", u"",
          u"1 proceso con p hilos, variante B (`--threads-per-core=1`%s), NELM = 4, en una sola cadena."
          % (u" + `srun --threads-per-core=1`" if geo.get("paso_tpc") == "1" else u""), u"",
          u"| Entorno | Variables | p | Mascara (nucleos) | Hilos activos (mediana) | Uso de CPU | Iteraciones | t/iter (s) | S por iteracion | Veredicto |",
          u"|---|---|---|---|---|---|---|---|---|---|"]
    for f in filas:
        md.append(u"| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s: %s |" % (
            f["E"], ENTORNO.get(f["E"], ("?",))[0], f["p"], f.get("nuc_mascara") if f.get("nuc_mascara") is not None else "-",
            "-" if f.get("activos") is None else "%g" % f["activos"], "-" if f.get("uso") is None else "%.2f" % f["uso"],
            "-" if f.get("n_iter") is None else f["n_iter"], "-" if not f.get("t_iter_s") else "%.1f" % f["t_iter_s"],
            "-" if not f.get("S_iter") else "%.2f" % f["S_iter"], f["veredicto"], f["motivo"]))
    if t_ref:
        md += [u"", u"Referencia con 1 hilo (A.1, p = 1, primeras 4 iteraciones): %.1f s por iteracion." % t_ref]
    ob = [f for f in filas if f.get("openblas_omp")]
    if ob:
        md += [u"OpenBLAS: %s; con OpenMP (enlaza libgomp): %s." % (ob[0].get("openblas_ver") or "?", ob[0]["openblas_omp"])]
    for f in filas:
        if f.get("pila"):
            md += [u"Pila de los hilos (gdb) en %s p=%d: %s." % (f["E"], f["p"], f["pila"])]
    if eleccion:
        md += [u"", u"**Entorno elegido para el barrido: %s** (%s)." % (eleccion, ENTORNO[eleccion][0])]
        if proy:
            md += [u"Amdahl con S(4) y S(24) por iteracion: f = %.3f; costo proyectado del barrido: %.0f h-nucleo."
                   % (proy["f"], proy["costo"])]
    with io.open(os.path.join(dd, "parteA2_hilos.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md) + u"\n")
    if eleccion and proy:
        lineas = [u"variante=B", u"paso_tpc=%s" % geo.get("paso_tpc", "0"), u"entorno=%s" % eleccion,
                  u"blas_hilos=%s" % ENTORNO[eleccion][1], u"espera_omp=%s" % ENTORNO[eleccion][2],
                  u"p_max=24", u"f_amdahl=%.4f" % proy["f"], u"t1_s=%.1f" % t1, u"costo_estimado_h=%.1f" % proy["costo"]]
        lineas += [u"limite_p%d=%s" % (p, proy["lim"][p]) for p in P_BARRIDO]
        lineas += [u"motivo=%s pasa con p = 4 y 24 sin sobresuscripcion" % eleccion]
        with io.open(CONFIG, "w", encoding="utf-8") as fh:
            fh.write(u"\n".join(lineas) + u"\n")
    elif os.path.exists(CONFIG):
        os.remove(CONFIG)


def guardar_evidencia(filas):
    base = os.path.join(REPO, "datos", "parteA2_hilos")
    for f in filas:
        d = f.get("dir")
        if not d or not os.path.isdir(d):
            continue
        dst = os.path.join(base, os.path.basename(d))
        os.path.isdir(dst) or os.makedirs(dst)
        for a in ("tiempo.txt", "entorno.txt", "afinidad.txt", "hilos.txt", "hilos_crudo.txt", "progreso.txt",
                  "slurm_job.txt", "carga.txt", "variables.txt", "topologia.txt", "pila_hilos.txt", "OSZICAR"):
            if os.path.isfile(os.path.join(d, a)):
                shutil.copy(os.path.join(d, a), os.path.join(dst, a))
        err = leer(os.path.join(d, "vasp.err")).splitlines()
        if err:
            with io.open(os.path.join(dst, "vasp_err_inicio.txt"), "w", encoding="utf-8") as fh:
                fh.write(u"\n".join(err[:60]) + u"\n")
        lp = [l for l in leer(os.path.join(d, "OUTCAR")).splitlines() if "LOOP:" in l or "running" in l]
        if lp:
            with io.open(os.path.join(dst, "outcar_loop.txt"), "w", encoding="utf-8") as fh:
                fh.write(u"\n".join(lp) + u"\n")
    out = git("add", "datos")
    if out.startswith("ERROR"):
        print("AVISO: git add fallo: " + out)
        return
    if git("diff", "--cached", "--quiet").startswith("ERROR"):
        out = git("commit", "-q", "-m", "Parte A.2: diagnostico de OpenMP, ronda 2 (OpenBLAS y espera) y entorno elegido")
        print(("AVISO: git commit fallo: " + out) if out.startswith("ERROR") else ("Commit hecho: " + git("log", "--oneline", "-1")))
    else:
        print("Nada nuevo que guardar en Git.")


def main():
    ap = argparse.ArgumentParser(description="Diagnostico de OpenMP, ronda 2 (A.2)")
    ap.add_argument("--guardar", action="store_true", help="guarda la evidencia en datos/ y hace commit")
    a = ap.parse_args()
    if not os.path.isfile(REG):
        print("Todavia no se envio la ronda 2: primero  bash slurm/parteA2_hilos.sh")
        return
    geo = clave_valor(GEO)
    t_ref, t1, n1 = referencia_a1()
    filas = evaluar(t_ref)
    pend = [f for f in filas if f["veredicto"] == "EN CURSO"]
    eleccion = None if pend else elegir(filas)
    proy = proyectar(filas, eleccion, t_ref, t1) if eleccion else None
    imprimir(filas, t_ref, t1, n1, eleccion, proy, geo)
    if a.guardar:
        print()
        if pend:
            print("No se guarda nada mientras haya corridas en cola o corriendo.")
            return
        guardar(filas, t_ref, t1, eleccion, proy, geo)
        guardar_evidencia(filas)
    elif not pend:
        guardar(filas, t_ref, t1, eleccion, proy, geo)


if __name__ == "__main__":
    main()
