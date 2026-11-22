# Examen final: benchmarking de paralelismo con VASP (Ruta A)

**Autor:** Rafael Silvera Altamiranda, Grupo de Física Aplicada, Universidad del Norte
**Curso:** Introducción a la Computación Científica
**Clúster:** Granado-HPC (Slurm), particiones `phys_fast` y `phys_general`
**Software:** VASP 6.3.0, módulo `vasp/6.3.0-gnu9-openmp`, binario `vasp_std`

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
   `AllocCPUS = p` (verificado en la A.1). Para hilos OpenMP, en cambio, `--cpus-per-task=p` entregó
   solo p/2 núcleos físicos (2026-11-18): la forma correcta se elige con el diagnóstico de la A.2.
   Cada script de envío pasa sus opciones de `--hint` o `--threads-per-core` de forma explícita.
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
   núcleo físico se comparta (`afinidad.txt`). Con OpenMP, además, mira con `ps -L` cada 10 s en qué
   núcleos físicos corren de verdad los hilos de VASP (`hilos.txt`). Si algo no se cumple, la corrida
   se marca como fallida.

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
    parteA2_diagnostico.sh          Tarea A.2: 6 corridas cortas para elegir como pedir los nucleos
    parteA2_enviar.sh               Tarea A.2: envia el barrido OpenMP en una cadena (o --faltantes)
    parteA2_estado.sh               Tarea A.2: estado del diagnostico y del barrido
  analisis/
    resumen_t0.py                   tiempos, memoria y verificaciones; elige y fija el caso base
    parteA1_resumen.py              Tarea A.1: t(p) con media y desviacion estandar, S(p), E(p)
    parteA2_diagnostico.py          Tarea A.2: veredicto del diagnostico de OpenMP y variante elegida
    parteA2_resumen.py              Tarea A.2: t(p) con media y desviacion estandar, S(p), E(p)
  datos/                            CSV, tabla del caso base y evidencia de la corrida elegida
  teoria/                           Parte I
  reporte/                          fuentes .tex y PDF
  (fuera de Git: POTCAR, lib/, corridas/, logs/)
```

Cada carpeta `corridas/<job>_<id>/` guarda las entradas y sus md5 (`entradas.md5`), `entorno.txt`
(nodo, procesador, commit), `modulos.txt`, `ldd.txt`, `variables.txt`, `slurm_job.txt`, la afinidad
(`vasp.err` y `afinidad.txt`), `carga.txt`, `tiempo.txt` y las salidas de VASP.

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

**Paso 1, diagnóstico.** Corridas cortas del caso base (`NELM = 4`) para elegir cómo pedir p núcleos
físicos para p hilos:

| Variante | Petición a Slurm | Probada con |
|---|---|---|
| A | `--ntasks=1 --cpus-per-task=p --hint=nomultithread` (la forma directa; falló el 2026-11-18) | p = 4 |
| B | `--ntasks=1 --cpus-per-task=p --threads-per-core=1` | p = 4 y 24 |
| C | `--ntasks=p --cpus-per-task=1 --hint=nomultithread` (la asignación de la A.1) y un solo proceso con `srun --ntasks=1 --cpus-per-task=p` | p = 4 y 24 |
| D | `--ntasks=1 --cpus-per-task=2p --hint=multithread` (2 hilos de hardware por núcleo: cobra el doble y solo cabe hasta p = 12) | p = 4 |

Una variante pasa si VASP corre con p hilos, la máscara del proceso tiene p núcleos físicos y los hilos
corren de verdad en p núcleos físicos distintos. Se prefiere la que además cobra p CPUs y llega a
p = 24 (C, luego B). El resultado queda en `datos/parteA2_diagnostico.md` y la variante elegida en
`datos/parteA2_variante.txt`.

**Paso 2, barrido.**

| Decisión | Valor | Justificación |
|---|---|---|
| Caso | `entradas/caso_base/` (K7), igual que en la A.1 | Solo cambia el número de hilos |
| Procesos e hilos | 1 proceso MPI con p hilos OpenMP, p = 2, 4, 8, 12, 16, 24; `OMP_PLACES=cores`, `OMP_PROC_BIND=close`, `OMP_STACKSIZE=512m` | Memoria compartida pura; variables recomendadas para el binario OpenMP de VASP |
| t(1) | Las 3 corridas con p = 1 de la A.1 | 1 proceso con 1 hilo es exactamente el mismo cálculo |
| Repeticiones | 3 por cada p (18 corridas), en una sola cadena | Igual que la A.1: sin competencia entre jobs propios |
| Orden | r1: 2→24; r2: 24→2; r3: 8, 16, 2, 24, 4, 12 | Que la hora no se confunda con el efecto de p |
| Verificaciones | Afinidad, hilos en p núcleos físicos, VASP con 1 proceso y p hilos, mismo trabajo, carga ajena, atípicas | Rúbrica: corrección técnica |

```bash
bash slurm/parteA2_preparar.sh                   # verifica y guarda los scripts en Git
bash slurm/parteA2_diagnostico.sh                # paso 1: 6 corridas cortas
bash slurm/parteA2_estado.sh                     # seguimiento
python3 analisis/parteA2_diagnostico.py          # veredicto y variante elegida
python3 analisis/parteA2_diagnostico.py --guardar
bash slurm/parteA2_enviar.sh                     # paso 2: barrido en una cadena
python3 analisis/parteA2_resumen.py              # t(p) media +/- desviacion, S(p), E(p)
python3 analisis/parteA2_resumen.py --guardar
```

---

## 7. Historial

El primer caso base (monocapa 3×3 con V_Se, ENCUT = 500 eV, calibrado el 2026-11-17) y las pruebas
del entorno de Granado quedaron en la carpeta `examen_silvera_monocapa/`, con su propio historial de
Git. Se reemplazó por el bulk 2×2×1 siguiendo la indicación del docente del 2026-11-21. Esas
corridas siguen contando en `sacct` y aparecen en `slurm/consumo.sh`.
