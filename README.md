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
~8 800 ondas planas por k-point, un SCF con la malla mínima (3×3×1, 3 k-points) dura solo unos
minutos. Por eso se calibra **solo la malla** (todo lo demás queda fijo) con siete candidatos que
corren a la vez, cada uno en 1 núcleo físico:

| Candidato | Malla Γ (supercelda) | Equivale en la primitiva | k-points irreducibles | t₁ estimado (min) |
|---|---|---|---|---|
| K1 | 3×3×1 | 6×6×1 | 3 | ~4 |
| K2 | 4×4×2 | 8×8×2 | 8 | ~11 |
| K3 | 6×6×2 | 12×12×2 | 14 | ~18 |
| K4 | 6×6×3 | 12×12×3 (densidad convergida en la tesis) | 17 | ~22 |
| K5 | 8×8×2 | 16×16×2 | 20 | ~26 |
| K6 | 8×8×3 | 16×16×3 | 25 | ~33 |
| K7 | 10×10×3 | 20×20×3 | 36 | ~47 |

Los k-points irreducibles se calcularon con la simetría del POSCAR (C₃ᵥ más inversión temporal) y
VASP debe reportar los mismos; `analisis/resumen_t0.py` lo verifica. Los tiempos son una estimación
a partir de la calibración real del intento anterior (monocapa 3×3: 45–88 s por iteración con
4 k-points y ~43 000 ondas planas) y **pueden variar un factor 2** según la carga del nodo.

**Regla de elección:** entre los candidatos que convergen con 15 ≤ t₁ ≤ 60 min, el más cercano a
**30 min** (centro de la ventana en escala logarítmica). Así queda un margen de un factor 2 hacia ambos
lados, que es la variabilidad que se midió en Granado entre corridas idénticas (20 vs 41 min).

Los candidatos también sirven para la **Tarea C.2**: dan el tiempo y la energía total en función de la
malla de k-points.

---

## 2. Configuración de Granado (verificada el 2026-11-17 y 18)

Estos hallazgos afectan a cualquier cálculo con VASP en Granado. `slurm/vasp_job.sh` ya los incorpora:

1. **Hyperthreading.** Slurm cuenta cada hilo como una CPU (`SelectTypeParameters = CR_CPU_MEMORY`).
   Con `--hint=nomultithread` cada proceso recibe un núcleo físico completo y `AllocCPUS = p`.
2. **Lanzador MPI.** `srun` sin opciones falla con OpenMPI 4.1.2 (`OPAL ERROR … pmix3x`) y `mpirun`
   ignora la asignación de Slurm (ata todos los procesos al núcleo 0). Se usa
   `srun --mpi=pmix --cpu-bind=verbose,cores`. Los avisos `PMIX ERROR: BAD-PARAM` son inofensivos;
   **no** se fija `PMIX_MCA_gds`, porque rompe `MPI_Init`.
3. **Librería FFTW.** Los nodos de cómputo no tienen `libfftw3_omp.so.3`. Se copia, junto con
   `libfftw3.so.3`, desde `/lib64` del login a `lib/`, y cada job la agrega a `LD_LIBRARY_PATH`.
4. **`/usr/bin/time` no existe en los nodos.** El tiempo de pared se mide con `date +%s.%N` alrededor
   de `srun`; también se guardan el `Elapsed time` del OUTCAR y el `Elapsed` de `sacct`.
5. **Memoria.** Sin `--mem`, Slurm pide la memoria de todo el nodo y el job queda bloqueado por
   `QOSMaxMemoryPerJob`. Se pide `--mem=4G` (el caso usa menos de 2 GB).
6. **Nodo `granado-compute-18`.** El 2026-11-18 mataba los jobs en 1–2 s; se excluye.
7. **Carga del socket.** Los jobs de otros usuarios en el mismo socket pueden duplicar el tiempo.
   Cada corrida registra cuántos CPUs del socket estaban ocupados (`carga.txt`).

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
  analisis/
    resumen_t0.py                   tiempos, memoria y verificaciones; elige y fija el caso base
  datos/                            CSV, tabla del caso base y evidencia de la corrida elegida
  teoria/                           Parte I
  reporte/                          fuentes .tex y PDF
  (fuera de Git: POTCAR, lib/, corridas/, logs/)
```

Cada carpeta `corridas/<job>_<id>/` guarda las entradas y sus md5 (`entradas.md5`), `entorno.txt`
(nodo, commit), `modulos.txt`, `ldd.txt`, `variables.txt`, `slurm_job.txt`, la afinidad (`vasp.err`),
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

## 5. Historial

El primer caso base (monocapa 3×3 con V_Se, ENCUT = 500 eV, calibrado el 2026-11-17) y las pruebas
del entorno de Granado quedaron en la carpeta `examen_silvera_monocapa/`, con su propio historial de
Git. Se reemplazó por el bulk 2×2×1 siguiendo la indicación del docente del 2026-11-21. Esas
corridas siguen contando en `sacct` y aparecen en `slurm/consumo.sh`.
