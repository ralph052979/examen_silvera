#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parteA2_hilos.py - Tarea A.2, ronda 2 del diagnostico de OpenMP: hilos de OpenBLAS y
                   politica de espera de OpenMP (diseno 2 x 2 con p = 4, y p = 24).

Para cada corrida de logs/jobs_A2h2.txt (intento 2; el 1 se resume) informa: forma de
lanzar el paso de srun, hilos ACTIVOS de VASP (>= 50 % de un CPU) frente a los nucleos
que tiene, uso total de CPU, tiempo por iteracion SCF y la aceleracion por iteracion
frente a 1 hilo (primeras 4 iteraciones de la A.1, p = 1). Concluye si E1 sirve para el
barrido (mascara correcta y sin hilos de mas con p = 4 y 24).

Uso (desde la carpeta del repositorio):
  python3 analisis/parteA2_hilos.py             # tabla, explicacion y conclusion
  python3 analisis/parteA2_hilos.py --guardar   # ademas guarda la evidencia en datos/ y hace commit

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
REG1 = os.path.join(REPO, "logs", "jobs_A2h.txt")      # intento 1 (variante B)
REG2 = os.path.join(REPO, "logs", "jobs_A2h2.txt")     # intento 2 (nucleos como la A.1 + paso automatico)
REG = REG2 if os.path.isfile(REG2) else REG1
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


def evaluar(t_ref, reg=None):
    filas = []
    for l in leer(reg or REG).splitlines():
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
        f["paso"] = tie.get("paso_variante") or ent.get("paso_variante") or "-"
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


def conclusion(filas):
    """E1 queda como entorno del barrido si sus corridas con p = 4 y 24 terminaron OK (mascara de p
    nucleos y sin hilos de mas), aunque sean lentas: medir esa lentitud es justamente la A.2."""
    a, b = ver(filas, "E1", 4), ver(filas, "E1", 24)
    ok = bool(a and b and a.get("estado") == "OK" and b.get("estado") == "OK")
    return {"ok": ok, "S4": a.get("S_iter") if a else None, "S24": b.get("S_iter") if b else None}


def imprimir(filas, t_ref, t1, n1, conc, intento1):
    print("Diagnostico de OpenMP, ronda 2 (Tarea A.2): hilos de OpenBLAS y espera de OpenMP")
    print()
    if intento1:
        print("Intento 1 (variante B: --threads-per-core=1): " + "; ".join(
            "%s p=%d %s" % (f["E"], f["p"], resumen_intento1(f)) for f in intento1) + ".")
        print("Intento 2 (este): nucleos pedidos como en la A.1 (--ntasks=p --hint=nomultithread) y paso de srun")
        print("elegido en cada corrida por su mascara real.")
        print()
    for E in ("E0", "E1", "E2", "E3"):
        print("  %s : %s" % (E, ENTORNO[E][0]))
    print()
    cab = "%-4s %-3s %-8s %-5s %-8s %-8s %-6s %-5s %-9s %-7s %-9s %s" % (
        "Ent", "p", "Job", "Paso", "Mascara", "Activos", "UsoCPU", "Iter", "t/iter", "S_iter", "Veredicto", "Detalle")
    print(cab)
    print("-" * len(cab))
    for f in filas:
        if f.get("t_iter_s"):
            ti = "%.1f s" % f["t_iter_s"]
        elif f.get("sin_iter_en"):
            ti = "> %d s" % f["sin_iter_en"]
        else:
            ti = "-"
        print("%-4s %-3s %-8s %-5s %-8s %-8s %-6s %-5s %-9s %-7s %-9s %s" % (
            f["E"], f["p"], f["jobid"] if f["jobid"] != "RECHAZADO" else "-", (f.get("paso") or "-")[:5],
            "-" if f.get("nuc_mascara") is None else "%d/%d" % (f["nuc_mascara"], f["p"]),
            "-" if f.get("activos") is None else "%g/%d" % (f["activos"], f["p"]),
            "-" if f.get("uso") is None else "%.2f" % f["uso"],
            "-" if f.get("n_iter") is None else f["n_iter"], ti,
            "-" if f.get("S_iter") is None else "%.2f" % f["S_iter"], f["veredicto"], f["motivo"]))
    print()
    print("Paso = forma de srun elegida por su mascara real (hint = --hint=nomultithread, tpc = --threads-per-core=1,")
    print("doble = 2 hilos de hardware por nucleo). Mascara = nucleos fisicos del proceso / hilos.")
    print("Activos = mediana de hilos con >= 50 % de un CPU.")
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
    if conc["ok"]:
        print("CONCLUSION: E1 (OPENBLAS_NUM_THREADS=1, espera por defecto) es correcto con p = 4 y 24: mascara de p")
        print("nucleos y sin hilos de mas. S por iteracion: S(4) = %s, S(24) = %s: OpenMP casi no acelera este caso."
              % ("-" if conc["S4"] is None else "%.2f" % conc["S4"], "-" if conc["S24"] is None else "%.2f" % conc["S24"]))
        print("El barrido usa E1 con una carga fija de 4 iteraciones SCF (el SCF completo costaria ~110 h-nucleo).")
        if os.path.isfile(os.path.join(REPO, "logs", "jobs_A2.txt")):
            print("(El barrido de la A.2 ya se envio y este diagnostico ya no pide ningun paso.)")
        else:
            print("Siguiente paso:  python3 analisis/parteA2_hilos.py --guardar   y despues   bash slurm/parteA2_enviar.sh")
    else:
        geo = [f for f in filas if (f.get("estado") or "").startswith("FALLA_geometria") and f.get("dir")]
        if geo:
            print("Pruebas del paso de srun en %s p=%d (geometria_pasos.txt):" % (geo[0]["E"], geo[0]["p"]))
            for l in leer(os.path.join(geo[0]["dir"], "geometria_pasos.txt")).splitlines():
                if "PMIX ERROR" not in l:
                    print("  " + l[:150])
        print("E1 no termino bien con p = 4 y p = 24: no envies el barrido y mandame una captura.")


def resumen_intento1(f):
    if f["veredicto"] == "OMITIDA":
        return "omitida"
    if f.get("nuc_mascara"):
        return "mascara %d/%d, %s" % (f["nuc_mascara"], f["p"], (f.get("estado") or "-").split("_(")[0].replace("FALLA_", ""))
    return f["veredicto"].lower()


def guardar(filas, t_ref, t1, conc, intento1):
    dd = os.path.join(REPO, "datos")
    os.path.isdir(dd) or os.makedirs(dd)
    cols = ["E", "p", "jobid", "nodo", "alloc", "nuc_mascara", "hilos_vasp", "hilos_totales", "activos", "uso",
            "n_iter", "t_iter_s", "cpu_real", "S_iter", "carga_ajena", "openblas_omp", "estado", "veredicto", "motivo"]
    with io.open(os.path.join(dd, "parteA2_hilos.csv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for f in filas:
            w.writerow(["" if f.get(c) is None else (("%.3f" % f[c]) if isinstance(f[c], float) else f[c]) for c in cols])
    md = [u"# Tarea A.2 - diagnostico de OpenMP, ronda 2 (hilos de OpenBLAS y espera de OpenMP)", u""]
    if intento1:
        md += [u"Intento 1, variante B (`--ntasks=1 --cpus-per-task=p --threads-per-core=1`): " + u"; ".join(
            u"%s p=%d %s" % (f["E"], f["p"], resumen_intento1(f)) for f in intento1) +
            u". La mascara se reviso al arrancar y las corridas se detuvieron en segundos.", u""]
    md += [u"Intento 2: 1 proceso con p hilos; nucleos pedidos como en la A.1 (`--ntasks=p --cpus-per-task=1 "
           u"--hint=nomultithread`); antes de VASP se prueban varias formas de lanzar el paso con srun y se usa la "
           u"primera cuya mascara real tiene p nucleos fisicos (`geometria_pasos.txt` de cada corrida). NELM = 4, "
           u"en una sola cadena.", u"",
          u"| Entorno | Variables | p | Paso de srun | Mascara (nucleos) | Hilos activos (mediana) | Uso de CPU | Iteraciones | t/iter (s) | S por iteracion | Veredicto |",
          u"|---|---|---|---|---|---|---|---|---|---|---|"]
    for f in filas:
        md.append(u"| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s: %s |" % (
            f["E"], ENTORNO.get(f["E"], ("?",))[0], f["p"], f.get("paso") or "-",
            f.get("nuc_mascara") if f.get("nuc_mascara") is not None else "-",
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
    if conc["ok"]:
        md += [u"", u"**Conclusion.** E1 (`OPENBLAS_NUM_THREADS=1`, espera de libgomp por defecto) es correcto con "
                    u"p = 4 y 24 (mascara de p nucleos fisicos, p hilos activos), pero OpenMP casi no acelera este caso: "
                    u"S por iteracion = %s con 4 hilos y %s con 24. La espera pasiva es mas lenta. El barrido de la A.2 "
                    u"usa E1 con una carga fija de 4 iteraciones SCF (README, seccion 6)."
               % ("-" if conc["S4"] is None else "%.2f" % conc["S4"], "-" if conc["S24"] is None else "%.2f" % conc["S24"])]
    with io.open(os.path.join(dd, "parteA2_hilos.md"), "w", encoding="utf-8") as fh:
        fh.write(u"\n".join(md) + u"\n")


def guardar_evidencia(filas):
    base = os.path.join(REPO, "datos", "parteA2_hilos")
    for f in filas:
        d = f.get("dir")
        if not d or not os.path.isdir(d):
            continue
        dst = os.path.join(base, os.path.basename(d))
        os.path.isdir(dst) or os.makedirs(dst)
        for a in ("tiempo.txt", "entorno.txt", "afinidad.txt", "hilos.txt", "hilos_crudo.txt", "progreso.txt",
                  "slurm_job.txt", "carga.txt", "variables.txt", "topologia.txt", "pila_hilos.txt",
                  "geometria_pasos.txt", ".afinidad_inicio.txt", "OSZICAR"):
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
        out = git("commit", "-q", "-m", "Parte A.2: diagnostico de OpenMP, ronda 2 (hilos de OpenBLAS y espera de OpenMP)")
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
    t_ref, t1, n1 = referencia_a1()
    filas = evaluar(t_ref)
    intento1 = evaluar(t_ref, REG1) if (REG == REG2 and os.path.isfile(REG1)) else []
    pend = [f for f in filas if f["veredicto"] == "EN CURSO"]
    conc = conclusion(filas) if not pend else {"ok": False, "S4": None, "S24": None}
    imprimir(filas, t_ref, t1, n1, conc, intento1)
    if a.guardar:
        print()
        if pend:
            print("No se guarda nada mientras haya corridas en cola o corriendo.")
            return
        guardar(filas, t_ref, t1, conc, intento1)
        guardar_evidencia(filas + intento1)
    elif not pend:
        guardar(filas, t_ref, t1, conc, intento1)


if __name__ == "__main__":
    main()
