# Examen final: benchmarking de paralelismo con VASP (Ruta A)

**Autor:** Rafael Silvera Altamiranda, Grupo de Física Aplicada, Universidad del Norte
**Curso:** Introducción a la Computación Científica
**Clúster:** Granado-HPC (Slurm), particiones `phys_fast` y `phys_general`
**Software:** VASP 6.3.0, módulo `vasp/6.3.0-gnu9-openmp`, binario `vasp_std`

**Reporte:** [`reporte/reporte.pdf`](reporte/reporte.pdf) (fuente `reporte/reporte.tex`). **Parte I:**
[`teoria/parte1.md`](teoria/parte1.md).

### Resultados principales

| Tarea | Resultado |
|---|---|
| Caso base | WSe₂ bulk 2×2×1 con V_Se (23 átomos), ENCUT = 400 eV, 10×10×3 (36 k-points): t(1) = 31.2 min |
| A.1 MPI | S(24) = 6.49 ± 0.13, E = 0.27; Amdahl f = 0.095, pero la fracción de Karp–Flatt crece de 0.08 a 0.12 |
| A.2 OpenMP | S máxima 1.39 con 4 hilos; con 24 hilos, 0.92 (binario sin `fftlib`) |
| B.1 débil | KPAR = p, 3 k-points por proceso: E_w(24) = 0.50, toda en el eigensolver (contención del nodo); memoria total ∝ p², pasa la QOS en p ≈ 29 |
| B.2 híbrido | óptimo 24 × 1 (89 s); 12 × 2 es 3 % más lento con 40 % menos memoria |
| C.1 | energías convergidas con p = 1…24: diferencias ≤ 2.6E-8 eV (7E5 ε_mach·\|E\|), 40 veces menos que EDIFF |
| C.2 | ENCUT = 400 eV es el punto de rendimientos decrecientes (0.33 meV/átomo) |
| D.1 | `phys_fast`, MPI puro, KPAR = 12 en 12 núcleos: 3.9 min y 0.79 h-núcleo por cálculo |
| Consumo | 148.5 de 250 h-núcleo (`datos/parteD_tabla.md`, `datos/sacct_registro_examen.txt`) |

---

## 1. Sistema físico y caso base (Tarea 0)

Indicación del docente (2026-11-21): usar el material de la tesis, **WSe₂ con vacancias**, en una
**supercelda 2×2×1 del bulk a la que se le quita un átomo de Se**, sin relajaciones (solo un SCF),
con una energía de corte no muy grande y una malla de k-points pequeña.

| Decisión | Valor | Justificación |
|---|---|---|
| Sistema | 2H-WSe₂ bulk, supercelda 2×2×1, una vacancia de Se (V_Se): **23 átomos** (8 W + 15 Se) | La V_Se es el defecto central de la tesis (el de menor energía de formación en todas las condiciones de potencial químico estudiadas) |
| Geometría | Experimental (Schutte, de Boer y Jellinek, *J. Solid State Chem.* 70, 207 (1987)): a = 3.282 Å, c = 12.96 Å, z(Se) = 0.621. W–Se = 2.527 Å. **No se relaja** | Indicación del docente. Además, todas las corridas hacen exactamente el mismo trabajo |
| Vacancia | Se en (1/6, 1/3, 0.621) de la supercelda. Quedan 3 W con 5 vecinos Se; simetría C₃ᵥ | En el bulk prístino todos los Se son equivalentes |
| PAW / funcional | `PAW_PBE W_sv` (14 e⁻) y `PAW_PBE Se` (6 e⁻), PBE; los mismos de la tesis. **NELECT = 202** | El POTCAR no va a Git (licencia). Su identidad queda en `entradas/POTCAR.titel` y `.sha256` |
| Cálculo | SCF estático desde cero: `ISTART = 0`, `ICHARG = 2`, `IBRION = -1`, `NSW = 0` | Sin relajación, como pidió el docente |
| `ENCUT` | **400 eV**, `PREC = Normal`, `LREAL = .FALSE.` | 20 % por debajo de los 500 eV convergidos en la tesis y 1.8 veces el ENMAX del W_sv (223 eV). `LREAL = .FALSE.` es lo que VASP recomienda para celdas pequeñas |
| Electrónica | `EDIFF = 1E-6`, `ALGO = Normal`, `ISMEAR = 0`, `SIGMA = 0.05`, `ISPIN = 1` | La V_Se no es magnética (resultado de la tesis) |
| `NBANDS` | **144, fijo** | Es múltiplo de 48 = mcm(1, 2, 4, 8, 12, 16, 24). Si no se fija, VASP redondea NBANDS para repartir las bandas (121 en 1 proceso, 144 en 24), y el trabajo cambiaría con p |
| Paralelismo de referencia | `NCORE = 1`, `KPAR = 1` | Se varían explícitamente en las Partes A y B |
| E/S | `LWAVE = LCHARG = .FALSE.` | Escribir archivos grandes añade ruido de disco que no es paralelismo |
| Semilla aleatoria | No aplica | VASP inicia las funciones de onda con `INIWAV = 1` y no expone una semilla (`RANDOM_SEED` solo afecta la dinámica molecular). La reproducibilidad se mide en la Tarea C.1 |

### Calibración de la malla de k-points

El enunciado exige que el caso base tarde **entre 15 min y 1 h en 1 núcleo**. Con 23 átomos y
~8 800 ondas planas por k-point, un SCF con la malla pequeña que sugirió el docente (3×3×1) dura solo
unos minutos. Por eso se calibró **solo la malla** (todo lo demás fijo) con siete candidatos en
1 núcleo físico cada uno (2026-11-21, `granado-compute-14`):

| Candidato | Malla Γ (supercelda) | Equivale en la primitiva | k-points irreducibles | t₁ estimado (min) | t₁ medido (min) |
|---|---|---|---|---|---|
| K1 | 3×3×1 | 6×6×1 | 3 | ~4 | 3.3 |
| K2 | 4×4×2 | 8×8×2 | 8 | ~11 | 8.4 |
| K3 | 6×6×2 | 12×12×2 | 14 | ~18 | 14.2 |
| K4 | 6×6×3 | 12×12×3 (densidad convergida en la tesis) | 17 | ~22 | 17.1 |
| K5 | 8×8×2 | 16×16×2 | 20 | ~26 | 19.9 |
| K6 | 8×8×3 | 16×16×3 | 25 | ~33 | 24.0 |
| **K7** | **10×10×3** | 20×20×3 | **36** | ~47 | **32.4** |

**Regla de elección (fijada antes de ver los datos):** entre los candidatos que convergen con
15 ≤ t₁ ≤ 60 min, el más cercano a **30 min** (centro de la ventana en escala logarítmica), para tener
margen de un factor 2 hacia ambos lados (variabilidad medida en Granado entre corridas idénticas:
20 vs 41 min).

**Resultado: caso base K7** (`entradas/caso_base/`), t₁ = 32.4 min, 18 iteraciones SCF, 913 MB
estimados por VASP y 937 MB medidos. Tabla completa en `datos/tarea0_tabla_caso_base.md`.

Observaciones:
- VASP usó en las 7 mallas exactamente los k-points irreducibles calculados con la simetría del
  POSCAR (C₃ᵥ más inversión temporal).
- Todas convergen en 18 iteraciones y el costo es proporcional a los k-points: ~2.9 s por k-point
  por iteración, más ~1 s fijo. Las estimaciones previas fueron ~25 % más altas porque suponían
  25 iteraciones.
- Los 36 k-points de K7 se reparten exactamente con KPAR = 1, 2, 3, 4, 6 y 12 en 24 núcleos, lo
  que permite estudiar KPAR en la Parte B (con 17 o 25 k-points no habría reparto parejo).
- K1–K5 corrieron los cinco a la vez en el mismo nodo; K7 corrió su último tramo solo y salió ~10 %
  más rápido por k-point: indicio de competencia por ancho de banda de memoria entre jobs.
- Los candidatos también sirven para la **Tarea C.2** (tiempo y energía en función de la malla).

---

## 2. Configuración de Granado (verificada el 2026-11-17 y 18)

Estos hallazgos afectan a cualquier cálculo con VASP en Granado. `slurm/vasp_job.sh` ya los incorpora:

1. **Hyperthreading.** Slurm cuenta cada hilo como una CPU (`SelectTypeParameters = CR_CPU_MEMORY`).
   Con `--ntasks=p --hint=nomultithread` cada proceso MPI recibe un núcleo físico completo y
   `AllocCPUS = p` (verificado en la A.1). Para 1 proceso con p hilos, en cambio, el paso de `srun`
   suele tomar los 2 hilos de hardware de p/2 núcleos: `--cpus-per-task=p` dio p/2 núcleos
   (2026-11-18), `srun -n 1 -c 24` dentro del job de la A.1 dio 12 núcleos (C24, 2026-11-22) y la
   variante B (`--threads-per-core=1`) dio 2 núcleos para 4 hilos (intento 1 de la ronda 2). Por eso,
   en la A.2 el job pide los núcleos igual que la A.1 y `vasp_job.sh` (`PASO_MODO=auto`) prueba, con
   un comando de 1 s, cuatro formas de lanzar el paso (`--hint=nomultithread`, `--threads-per-core=1`,
   sin opciones y con 2p CPUs) y usa la primera cuya máscara **real** tiene p núcleos físicos
   (`geometria_pasos.txt`). Cada script de envío pasa sus opciones de forma explícita.
2. **Lanzador MPI.** `srun` sin opciones falla con OpenMPI 4.1.2 (`OPAL ERROR … pmix3x`) y `mpirun`
   ignora la asignación de Slurm (ata todos los procesos al núcleo 0). Se usa
   `srun --mpi=pmix --cpu-bind=verbose,cores`. Los avisos `PMIX ERROR: BAD-PARAM` son inofensivos;
   **no** se fija `PMIX_MCA_gds`, porque rompe `MPI_Init`.
3. **Librería FFTW.** Los nodos de cómputo no tienen `libfftw3_omp.so.3`. Se copia, junto con
   `libfftw3.so.3`, desde `/lib64` del login a `lib/`, y cada job la agrega a `LD_LIBRARY_PATH`.
4. **`/usr/bin/time` no existe en los nodos.** El tiempo de pared se mide con `date +%s.%N` alrededor
   de `srun`; también se guardan el `Elapsed time` del OUTCAR y el `Elapsed` de `sacct`.
5. **Memoria.** Sin `--mem`, Slurm pide la memoria de todo el nodo y el job queda bloqueado por
   `QOSMaxMemoryPerJob`. Se pide siempre de forma explícita: 4 GB en la Tarea 0 y 4 + p/2 GB en la A.1.
6. **Nodo `granado-compute-18`.** El 2026-11-18 mataba los jobs en 1–2 s; se excluye.
7. **Carga del socket.** Los jobs de otros usuarios en el mismo socket pueden duplicar el tiempo.
   Cada corrida registra cada minuto cuántos CPUs de cada socket están ocupados (`carga.txt`) y
   resume la carga ajena en los sockets del job (`carga_ajena_max` en `tiempo.txt`).
8. **Máximo 5 jobs corriendo por usuario** (`QOSMaxJobsPerUserLimit`, visto el 2026-11-21). Los
   demás esperan en cola aunque haya núcleos libres.
9. **Afinidad.** Cada corrida verifica, con las máscaras de `srun --cpu-bind=verbose` y la topología
   de `/sys`, que cada proceso tenga exactamente `--cpus-per-task` núcleos **físicos** y que ningún
   núcleo físico se comparta (`afinidad.txt`). La máscara se revisa al arrancar: si está mal, la corrida
   se detiene en segundos. Con OpenMP, además, cada 10 s se lee `/proc/<pid>/task/*/stat` de VASP (el
   pid lo escribe `srun` en `vasp.err`): en qué núcleos físicos están los hilos y cuántos están
   **activos** (≥ 50 % de un CPU) (`hilos_crudo.txt`, `hilos.txt`). Si hay más hilos activos que
   núcleos, la corrida se marca como sobresuscrita. (El primer monitor buscaba el proceso con `ps` y
   `/proc/<pid>/environ` y en la ronda 1 de la A.2 no obtuvo ninguna muestra.)
10. **OpenMP: un hilo de hardware por núcleo, o nada (A.2, 2026-11-22).** En la ronda 1, las corridas
    con p CPUs para p hilos no terminaron 4 iteraciones en 20 min; durante C24, `ps -L` mostró los 24 hilos
    de OpenMP en 12 núcleos físicos (2 por núcleo) y un hilo más trabajando. Con 4 núcleos pero con sus
    2 hilos de hardware en la máscara (D4), cada iteración tardó 127–133 s. En la ronda 2,
    con 1 hilo de hardware por núcleo, la misma iteración con 4 hilos tardó 69–70 s y trabajaron exactamente
    4 hilos. `OPENBLAS_NUM_THREADS=1` no cambia nada (E0 = 69,1 s, E1 = 69,8 s): la hipótesis de un hilo
    propio de OpenBLAS no se confirmó (la librería no enlaza libgomp). La espera pasiva
    (`OMP_WAIT_POLICY=PASSIVE`) es alrededor de 50 % más lenta con 4 hilos y 67 % con 24. Detalle en
    `datos/parteA2_diagnostico.md` y `datos/parteA2_hilos.md`.
11. **Salida de VASP con bufer.** VASP compilado con gfortran escribe OSZICAR y OUTCAR en bloques y los
    vacía al terminar: una corrida que Slurm corta los deja vacíos o incompletos. `vasp_job.sh` fija
    `GFORTRAN_UNBUFFERED_PRECONNECTED=y`, de modo que `vasp.out` muestra cada iteración SCF al momento,
    y anota en `progreso.txt` el segundo en que termina cada una. Además: si a los `MAX_PRIMERA_ITER`
    segundos no hay ninguna iteración, la corrida se detiene; y si Slurm avisa que se acaba el tiempo
    (`--signal=B:USR1@90`), se detiene VASP y se completan igual las verificaciones.

---

## 3. Estructura

```
examen_silvera/
  README.md                         este archivo
  entradas/
    generar_poscar_bulk.py          construye la supercelda 2x2x1 con V_Se y verifica la geometria
    POSCAR_bulk_2x2x1_VSe           supercelda generada (23 atomos)
    geometria_bulk_2x2x1_VSe.txt    a, c, z, distancias, Se removido
    INCAR_base                      INCAR del SCF (igual para todos los candidatos)
    candidatos.txt                  mallas de k-points de la calibracion
    candidatos/K1..K7/              INCAR, KPOINTS y POSCAR de cada candidato
    caso_base/                      caso elegido (lo crea analisis/resumen_t0.py --fijar)
    debil/P1..P24/                  Parte B.1: KPAR = p y 3p k-points (lo crea slurm/parteBC_preparar.sh)
    encut/E250..E800/               Tarea C.2: ENCUT = 250...800 eV con la malla 3x3x1 (idem)
    POTCAR.titel, POTCAR.sha256     identidad del POTCAR (el POTCAR no va a Git)
  slurm/
    preparar.sh                     prepara entradas, librerias y Git (login, sin VASP)
    enviar.sh                       envia los candidatos con sbatch
    vasp_job.sh                     job generico: una corrida de VASP en corridas/<job>_<id>/
    estado.sh                       estado de los jobs
    nuevo_candidato.sh              agrega y envia otra malla (solo si el resumen lo pide)
    consumo.sh                      horas-nucleo consumidas (sacct)
    parteA1_preparar.sh             Tarea A.1: verifica el caso base y guarda los scripts en Git
    parteA1_enviar.sh               Tarea A.1: envia las 21 corridas en una cadena (o --faltantes)
    parteA1_estado.sh               Tarea A.1: estado de las corridas
    parteA2_preparar.sh             Tarea A.2: verifica y guarda los scripts en Git
    parteA2_diagnostico.sh          Tarea A.2, ronda 1: 6 corridas cortas (como pedir los nucleos)
    parteA2_hilos.sh                Tarea A.2, ronda 2: hilos de OpenBLAS y espera de OpenMP
    parteA2_enviar.sh               Tarea A.2: barrido OpenMP con carga fija (4 iteraciones), en una cadena
    parteA2_estado.sh               Tarea A.2: estado de las dos rondas y del barrido
    parteBC_preparar.sh             Partes B y C.2: genera entradas/debil/ y entradas/encut/, verifica y guarda en Git
    parteB_enviar.sh                Parte B: B.1 (21 corridas) y B.2 (19) en una cadena (o --faltantes)
    parteC_enviar.sh                Tarea C.2: 3 series de ENCUT, un job de 1 nucleo por serie (o --faltantes)
    parteC2_serie.sh                Tarea C.2: job que corre vasp_job.sh una vez por cada ENCUT de la serie
    parteBC_estado.sh               Partes B y C.2: estado de las corridas
    empaquetar_datos.sh             junta logs/, datos/, sacct y los archivos pequenos de corridas/ (sin POTCAR)
    finalizar.sh                    analisis final: corre los scripts de analisis/ y guarda todo en Git
    publicar.sh                     agrega a Git los archivos pequenos de corridas/ y sube el repositorio
  analisis/
    resumen_t0.py                   tiempos, memoria y verificaciones; elige y fija el caso base
    parteA1_resumen.py              Tarea A.1: t(p) con media y desviacion estandar, S(p), E(p)
    parteA2_diagnostico.py          Tarea A.2, ronda 1: mascaras (tambien de las cortadas) y tiempos
    parteA2_hilos.py                Tarea A.2, ronda 2: hilos activos, tiempos por iteracion y conclusion
    parteA2_resumen.py              Tarea A.2: t(p), S(p), E(p) y comparacion por iteracion con MPI (A.1)
    parteA3_graficas.py             Tarea A.3: figuras S(p) y E(p), ajustes de Amdahl y de Karp-Flatt
    parteB_resumen.py               Parte B: Ew(p) del escalamiento debil; t, S, E y memoria del hibrido
    parteB_graficas.py              Parte B: figuras B1_debil y B2_hibrido, tiempo por rutina y modelo de memoria
    parteC_analisis.py              Parte C: energias con p (C.1), demostracion de la reduccion, ENCUT (C.2)
    parteD_colas.py                 Parte D: consumo por parte, espera en cola y costo por calculo (sacct)
    reporte_tablas.py               tablas del reporte (reporte/tablas/*.tex) a partir de datos/
  datos/                            CSV, tabla del caso base y evidencia de la corrida elegida
  teoria/                           Parte I: parte1.md (desarrollo), parte1_calculos.py y su salida
  reporte/                          reporte.tex, reporte.pdf, figuras/ y tablas/ (las generan los scripts de analisis/)
  corridas/<job>_<id>/              de cada corrida, solo los archivos pequenos (tiempos, afinidad, hilos,
                                    INCAR, KPOINTS, OSZICAR, OUTCAR): el analisis se repite desde un clon
  logs/jobs_*.txt                   listas de jobs de cada parte (que corrida es que job)
  (fuera de Git: POTCAR, lib/, archivos grandes de VASP y las salidas .out de Slurm)
```

Cada carpeta `corridas/<job>_<id>/` guarda las entradas y sus md5 (`entradas.md5`), `entorno.txt`
(nodo, procesador, commit, entorno de hilos, OpenBLAS), `modulos.txt`, `ldd.txt`, `variables.txt`,
`slurm_job.txt`, `topologia.txt`, la afinidad (`vasp.err` y `afinidad.txt`), los hilos (`hilos_crudo.txt`
y `hilos.txt`, que desde la Parte B incluye la memoria maxima VmHWM de cada proceso), `progreso.txt`,
`carga.txt`, `tiempo.txt` y las salidas de VASP.

---

## 4. Cómo reproducir la Tarea 0

```bash
cd /work/physics/users/rsilvera/examen_silvera
bash slurm/preparar.sh                       # login: POTCAR, POSCAR, candidatos, librerias, Git
bash slurm/enviar.sh                         # envia K1..K7 (1 nucleo cada uno)
bash slurm/estado.sh                         # seguimiento
python3 analisis/resumen_t0.py               # tabla, verificaciones y recomendacion
python3 analisis/resumen_t0.py --fijar KX    # fija el caso base, escribe la tabla y hace commit
bash slurm/consumo.sh                        # horas-nucleo usadas
```

Todo el cómputo se hace con `sbatch`. En el login solo se copian y verifican archivos, se genera el
POSCAR (menos de 1 s) y se leen resultados.

La tabla de la Tarea 0 queda en `datos/tarea0_tabla_caso_base.md` y la calibración completa en
`datos/tarea0_calibracion.csv`.

---

## 5. Parte A.1: escalamiento fuerte con MPI

| Decisión | Valor | Justificación |
|---|---|---|
| Caso | `entradas/caso_base/` (K7), el mismo INCAR, KPOINTS y POSCAR en todas las corridas | Solo cambia p |
| Procesos | p = 1, 2, 4, 8, 12, 16, 24 procesos MPI, 1 hilo cada uno, un núcleo físico por proceso (`--hint=nomultithread`) | Memoria distribuida pura |
| Parámetros de VASP | `NCORE = 1`, `KPAR = 1` (valores por defecto): VASP reparte las 144 bandas entre los p procesos (NPAR = p) | NCORE y KPAR se exploran en la Parte B; aquí se mide el efecto de p |
| Repeticiones | 3 por cada p (21 corridas) | Media y desviación estándar, como pide el enunciado |
| Ejecución | **Una sola cadena** (`--dependency=afterany`): nunca corren dos corridas propias a la vez | En la Tarea 0 se vio que los jobs propios en el mismo nodo se frenan entre sí (~10 %) |
| Orden | r1: 1→24; r2: 24→1; r3: 8, 1, 16, 4, 24, 2, 12 | Que la carga del clúster a distintas horas no se confunda con el efecto de p |
| Tiempo medido | Tiempo de pared de `srun` (`date` antes y después); también se guardan el `Elapsed` de VASP y el tiempo por iteración | El enunciado pide el tiempo de pared |
| Verificaciones | Afinidad por núcleo físico, mismo número de iteraciones, k-points y bandas en todas las corridas, carga ajena, corridas atípicas (> 15 % de la mediana) | Rúbrica: diseño experimental y corrección técnica |
| Memoria pedida | 4 + p/2 GB | En 1 proceso el caso usa 0.94 GB; cada proceso extra agrega ~0.3 GB |
| Límite de tiempo | 1:40 h (p = 1), 1 h (2), 40 min (4), 30 min (8), 20 min (12–24) | Unas 3 veces el tiempo esperado |

Las energías de todas las corridas quedan en `datos/parteA1_corridas.csv` con todos los dígitos que
imprime VASP: las de p = 1, 4, 16 y 24 sirven para la **Tarea C.1** sin correr nada más.

```bash
bash slurm/parteA1_preparar.sh            # verifica y guarda los scripts en Git
bash slurm/parteA1_enviar.sh              # envia la cadena de 21 corridas
bash slurm/parteA1_estado.sh              # seguimiento
python3 analisis/parteA1_resumen.py       # t(p) media +/- desviacion, S(p), E(p)
python3 analisis/parteA1_resumen.py --guardar   # evidencia a datos/ y commit
bash slurm/parteA1_enviar.sh --faltantes  # solo si alguna corrida fallo o se excluyo
```

Para excluir una corrida contaminada (por ejemplo, con carga ajena alta y tiempo atípico) se agrega
una línea `JOBID motivo` a `datos/parteA1_excluidas.txt`; `--faltantes` la repite.

---

## 6. Parte A.2: escalamiento fuerte con OpenMP

**Paso 1, diagnóstico ronda 1: cómo pedir p núcleos físicos.** Corridas cortas del caso base
(`NELM = 4`), enviadas a la vez:

| Variante | Petición a Slurm | Probada con |
|---|---|---|
| A | `--ntasks=1 --cpus-per-task=p --hint=nomultithread` | p = 4 |
| B | `--ntasks=1 --cpus-per-task=p --threads-per-core=1` | p = 4 y 24 |
| C | `--ntasks=p --cpus-per-task=1 --hint=nomultithread` (la asignación de la A.1) y un solo proceso con `srun --ntasks=1 --cpus-per-task=p` | p = 4 y 24 |
| D | `--ntasks=1 --cpus-per-task=2p --hint=multithread` (2 hilos de hardware por núcleo: cobra el doble) | p = 4 |

Resultado (hallazgos 1 y 10 de la sección 2): ninguna variante sirve tal cual. Tabla, máscaras leídas
de `vasp.err` y tiempos por iteración en `datos/parteA2_diagnostico.md`.

**Paso 2, diagnóstico ronda 2: hilos de OpenBLAS y espera de OpenMP.** El intento 1 usó la
variante B: la máscara de E1 con p = 4 tuvo 2 núcleos físicos y la corrida se detuvo al arrancar
(las demás, igual u omitidas). El intento 2 pidió los núcleos como en la A.1 y eligió el paso de `srun`
por su máscara real (hallazgo 1; en todas las corridas sirvió `--hint=nomultithread`). Diseño 2 × 2 con
p = 4 y, con p = 24, las que pasaron con p = 4; en una sola cadena:

| Entorno | Variables | t/iteración, p = 4 | t/iteración, p = 24 |
|---|---|---|---|
| E0 | las de la ronda 1 (control) | 69,1 s | — |
| E1 | `OPENBLAS_NUM_THREADS=1` | 69,8 s | 129,6 s |
| E2 | `OMP_WAIT_POLICY=PASSIVE` | 103,9 s | — |
| E3 | ambas | 103,2 s | 216,2 s (no terminó 4 iteraciones en 15 min) |

Con 1 hilo (A.1, p = 1), cada una de las primeras 4 iteraciones tarda 95,8 s: con OpenMP la iteración
se acelera 1,37 veces con 4 hilos y es **más lenta** con 24 (0,74). El binario de VASP no contiene la
cadena `fftlib` (`cadena_fftlib_en_vasp=0`), es decir, no se compiló con fftlib; su documentación advierte que sin fftlib el costo de planificar las FFT con
hilos se vuelve prohibitivo con muchos hilos, y que con OpenMP los hilos reemplazan a NCORE (se reparten
el trabajo de cada orbital; nuestra malla FFT es de 32×32×64). Tabla completa en `datos/parteA2_hilos.md`.

**Paso 3, barrido con carga fija.** Con esa eficiencia, cada corrida completa con p hilos cuesta casi
p × 30 min: el barrido completo (3 repeticiones, p = 2…24) costaría unas 110 h‑núcleo. Se mide en cambio
una carga fija de 4 iteraciones SCF del mismo caso:

| Decisión | Valor | Justificación |
|---|---|---|
| Caso | `entradas/caso_base/` (K7) con `NELM = 4` | Mismo cálculo que la A.1; carga fija e idéntica para todo p |
| Procesos e hilos | 1 proceso con p hilos, p = 1, 2, 4, 8, 12, 16, 24; núcleos pedidos como en la A.1 y paso de `srun` elegido por su máscara real; `OMP_PLACES=cores`, `OMP_PROC_BIND=close`, `OMP_STACKSIZE=512m`, `OPENBLAS_NUM_THREADS=1` | Memoria compartida pura; mismos núcleos y mismo cobro que la A.1; entorno del paso 2 |
| t(1) | Las 3 corridas con p = 1 de este barrido | Misma carga que las demás |
| Comparación con MPI | Tiempo por iteración (`LOOP` del OUTCAR): primeras 4 iteraciones de la A.1 frente a las 4 de la A.2 | Mismo trabajo por iteración en ambos modos |
| Repeticiones | 3 por cada p (21 corridas), en una sola cadena | Igual que la A.1 |
| Orden | r1: 1→24; r2: 24→1; r3: 8, 16, 1, 24, 4, 12, 2 | Que la hora no se confunda con el efecto de p |
| Costo | unas 31 h‑núcleo | Dentro de lo sugerido para la Parte A (90 h) |
| Verificaciones | Afinidad, hilos en p núcleos físicos, sin hilos activos de más, VASP con 1 proceso y p hilos, 4 iteraciones, carga ajena, atípicas | Rúbrica: corrección técnica |

```bash
bash slurm/parteA2_preparar.sh                   # verifica y guarda los scripts en Git
python3 analisis/parteA2_diagnostico.py --guardar   # ronda 1: evidencia
bash slurm/parteA2_hilos.sh                      # ronda 2 (intento 2): 6 corridas cortas en cadena
bash slurm/parteA2_estado.sh                     # seguimiento
python3 analisis/parteA2_hilos.py --guardar      # ronda 2: tabla, conclusion y evidencia
bash slurm/parteA2_enviar.sh                     # barrido con carga fija, en una cadena
python3 analisis/parteA2_resumen.py              # t(p), S(p), E(p) y comparacion por iteracion con MPI
python3 analisis/parteA2_resumen.py --guardar
```

---

## 7. Tarea A.3: gráficas y ajustes

`analisis/parteA3_graficas.py` lee las tablas de la A.1 y la A.2 y hace:

* `reporte/figuras/A3_escalamiento.pdf` (y `.png`): S(p) y E(p) de MPI y OpenMP con el tiempo de la
  corrida, la curva ideal, la ley de Amdahl ajustada (línea continua) y el modelo con sobrecosto
  (punteada).
* `reporte/figuras/A3_por_iteracion.pdf`: S(p) por iteración SCF de ambos modos y la fracción serial
  de Karp-Flatt e(p) = (1/S − 1/p)/(1 − 1/p) de cada punto.
* `datos/parteA3_ajustes.md` y `.csv`, `datos/parteA3_karp_flatt.csv`.

Modelos: Amdahl, 1/S = f + (1 − f)/p; con sobrecosto, 1/S = f + (1 − f)/p + c(p − 1). Los dos son
lineales en sus parámetros sobre y = 1/S − 1/p, y se ajustan por mínimos cuadrados ponderados con la
incertidumbre de S (p > 1). Si e(p) crece con p, el límite no es una fracción serial fija sino un
sobrecosto que crece con p (comunicación, sincronización, ancho de banda de memoria).

```bash
python3 analisis/parteA3_graficas.py --guardar   # necesita matplotlib (python3 -m pip install --user matplotlib)
```

---

## 8. Parte B: escalamiento débil e híbrido

| Decisión | B.1 (escalamiento débil) | B.2 (híbrido MPI × OpenMP) |
|---|---|---|
| Caso | `entradas/debil/P<p>/`: caso base con `KPAR = p` y una malla Gamma con 3p k-points irreducibles (3×3×1, 3×3×2, 5×5×3, 6×6×4, 10×10×3, 14×14×2, 10×10×6 para p = 1…24) | `entradas/caso_base/` (K7, 36 k-points) con `NELM = 4`, como la A.2 |
| Cómo crece el problema | Cada proceso calcula 3 k-points con las 144 bandas: el trabajo por proceso es fijo | — (24 núcleos fijos) |
| Procesos e hilos | p procesos MPI, p = 1, 2, 4, 8, 12, 16, 24; núcleos y `srun` como en la A.1 | 24 núcleos pedidos como en la A.1; paso de `srun` con N procesos × T hilos: (24,1) (12,2) (8,3) (6,4) (4,6) (2,12) (1,24), elegido por su máscara real (cada proceso con T núcleos físicos propios); `OMP_PLACES=cores`, `OMP_PROC_BIND=close` |
| Métrica | Ew(p) = t(1)/t(p), de la corrida y por iteración SCF (el número de iteraciones cambia con la malla) | t, S = t(1)/t con t(1) = p = 1 de la A.2 (misma carga), E = S/24 |
| Memoria | VmHWM de cada proceso cada 3 s: máximo por proceso y suma | Igual; para las corridas de la A.2 (1 proceso) se usa el OUTCAR |
| Repeticiones | 3 por p (21 corridas) | 3 por combinación (18) y 1 nueva de (1,24), que suma las 3 de la A.2 |
| Orden | Una sola cadena con B.1 y B.2 intercaladas por repetición, en orden distinto en cada una | |
| Entorno | `OPENBLAS_NUM_THREADS=1` en todas (en la A.2 no cambió el tiempo) | |
| Costo estimado | unas 17 h‑núcleo | unas 22 h‑núcleo |
| Verificaciones | Afinidad, SCF convergido, NKPTS = 3p, KPAR = p, p procesos | Afinidad (N procesos × T núcleos físicos, sin compartir), hilos en N×T núcleos, VASP con N procesos y T hilos, 4 iteraciones; procesos repartidos en los 2 sockets (B.3) |

```bash
bash slurm/parteBC_preparar.sh                   # genera entradas/debil/ y entradas/encut/, verifica y guarda en Git
bash slurm/parteB_enviar.sh                      # 40 corridas en una cadena
bash slurm/parteC_enviar.sh                      # C.2: 3 series de ENCUT, despues de la Parte B
bash slurm/parteBC_estado.sh                     # seguimiento
python3 analisis/parteB_resumen.py               # tablas B.1 y B.2, combinacion optima
python3 analisis/parteB_resumen.py --guardar     # ademas copia la evidencia a datos/parteB/ y hace commit
```

---

## 9. Tarea C.2: precisión frente a costo (ENCUT)

| Decisión | Valor | Justificación |
|---|---|---|
| Parámetro | ENCUT = 250, 300, 350, 400, 450, 500, 550, 600, 700, 800 eV (`entradas/encut/`) | Controla el error de truncamiento de la base de ondas planas |
| Referencia | ENCUT = 800 eV | Error de cada ENCUT: \|E(ENCUT) − E(800)\| por átomo |
| Malla | Gamma 3×3×1 (3 k-points) | El error de la base por átomo casi no depende de la malla; con 3 k-points la serie cuesta unas 3 h‑núcleo con 3 repeticiones (con la malla del caso base costaría unas 12 veces más) |
| Núcleos | 1 (como la p = 1 de la A.1), SCF completo | Tiempo de cómputo = horas‑núcleo |
| Repeticiones | 3 series, un job de 1 núcleo por serie, en orden distinto cada una (ascendente, descendente, mezclado) | Barras de error del tiempo; solo 3 jobs en la cola |
| Precisión simple | VASP 6.3.0 para CPU solo calcula en doble precisión | Se discute con ε_mach en la C.1 |

---

## 10. Análisis final y reporte (A.3, B, C, D)

Todo el análisis corre en el login con Python (sin VASP) a partir de `corridas/`, `logs/` y `datos/`:

```bash
bash slurm/finalizar.sh      # corre los scripts de abajo, verifica y guarda todo en Git
```

| Script | Escribe |
|---|---|
| `analisis/parteB_resumen.py` | `datos/parteB1_*.csv`, `parteB2_*.csv`, `parteB_tabla.md` |
| `analisis/parteA3_graficas.py` | `reporte/figuras/A3_*.pdf/png`, `datos/parteA3_*` |
| `analisis/parteB_graficas.py` | `reporte/figuras/B1_debil`, `B2_hibrido`, `datos/parteB1_rutinas.csv`, `parteB1_memoria.md` |
| `analisis/parteC_analisis.py` | `reporte/figuras/C2_precision_costo`, `datos/parteC1_*.csv`, `parteC2_*.csv`, `parteC_tabla.md` |
| `analisis/parteD_colas.py` | `datos/parteD_*.csv`, `parteD_tabla.md`; con `sacct`, `datos/sacct_colas.txt` y `datos/sacct_registro_examen.txt` (el comando del enunciado) |
| `analisis/reporte_tablas.py` | `reporte/tablas/*.tex` (ningún número de las tablas del reporte se copia a mano) |
| `teoria/parte1_calculos.py` | `teoria/parte1_salida.txt` |

El reporte se compila desde `reporte/` con `pdflatex reporte.tex` (dos veces) o en Overleaf subiendo
`reporte.tex`, `figuras/` y `tablas/`. Usa `lmodern`, `microtype`, `booktabs`, `adjustbox`, `caption`,
`enumitem` y `hyperref`; si `babel` tiene español, lo usa.

Decisiones del análisis:

* **C.1** usa la última línea `DAV` del OSZICAR (12 cifras significativas, resolución 1E-9 eV) y expresa
  las diferencias como |E(p) − E(1)| / (ε_mach |E(1)|), con ε_mach = 2⁻⁵². Compara también la energía
  después de exactamente 4 iteraciones de todas las corridas de la A.1, la A.2 y la B.2 (depende solo del
  número de procesos MPI) y muestra, con NumPy, el efecto del orden de una reducción en precisión simple y
  doble (semilla fija).
* **C.2** toma como referencia ENCUT = 800 eV y el error por átomo |E − E(800)| / 23; la serie usa la malla
  3×3×1 (el error de la base casi no depende de la malla) y 1 núcleo, con 3 repeticiones.
* **D** mide la espera en cola como `Start − Eligible` (en una cadena, `Eligible` es cuando termina el job
  anterior) y el costo por cálculo como p·t del SCF completo del caso base.

---

## 11. Historial

El primer caso base (monocapa 3×3 con V_Se, ENCUT = 500 eV, calibrado el 2026-11-17) y las pruebas
del entorno de Granado quedaron en la carpeta `examen_silvera_monocapa/`, con su propio historial de
Git. Se reemplazó por el bulk 2×2×1 siguiendo la indicación del docente del 2026-11-21. Esas
corridas siguen contando en `sacct` y aparecen en `slurm/consumo.sh`.
