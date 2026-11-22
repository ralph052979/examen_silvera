# Tarea 0 - caso base

| Propiedad | Valor |
|---|---|
| Sistema fisico | 2H-WSe2 bulk, supercelda 2x2x1 con una vacancia de Se (V_Se), el defecto de mi tesis |
| Geometria | experimental (Schutte et al. 1987): a = 3.2820 A, c = 12.9600 A, z(Se) = 0.6210; sin relajar |
| Numero de atomos | 23 (8 W + 15 Se) |
| Electrones de valencia (NELECT) | 202 |
| PAW / funcional | PAW_PBE W_sv 04Sep2015; PAW_PBE Se 06Sep2000 / PBE (GGA = PE) |
| Calculo | SCF estatico desde cero (ISTART = 0, ICHARG = 2, IBRION = -1, NSW = 0), EDIFF = 1E-6, ISMEAR = 0, SIGMA = 0.05, ISPIN = 1 |
| ENCUT / PREC / LREAL | 400 eV / Normal / .FALSE. |
| Malla de k-points | Gamma 10x10x3 (= 20x20x3 en la primitiva) -> 36 k-points irreducibles |
| NBANDS | 144 (fijo, multiplo de 48: no cambia con el numero de procesos) |
| Ondas planas por k-point (max) | 8815 |
| Malla FFT (NGX x NGY x NGZ) | 32 x 32 x 64 |
| Iteraciones SCF | 18 |
| t(1 nucleo) | 32.4 min (pared 1942 s; Elapsed de VASP 1942 s) |
| Memoria estimada por VASP | 913 MB |
| Memoria medida | 937 MB (OUTCAR) / 934 MB (sacct MaxRSS) |
| Limite de memoria de la QOS | 62 600 MB |
| Software | vasp.6.3.0 20Jan22 (build Nov 16 2026 16:16:35) complex; modulos gnu9/9.3.0, fftw/3.3.8-gnu9, openmpi4/4.1.2, wannier90/3.1.0-gnu9-openblas, openblas/0.3.17-gnu9, libevent/2.1.12, hdf5/1.10.6-gnu9, vasp/6.3.0-gnu9-openmp, scalapack/2.1.0; binario vasp_std |
| Nodo / CPUs ocupados en su socket | granado-compute-14.granado.uninorte.edu.co / 5/36 |
| Semilla aleatoria | No aplica: VASP inicia las funciones de onda con INIWAV = 1 y no expone una semilla (RANDOM_SEED solo afecta la dinamica molecular). La reproducibilidad se mide en la C.1 |
| Job / commit de las entradas | 529167 / 265ac0b |

Calibracion completa: `datos/tarea0_calibracion.csv`. Evidencia de la corrida: `datos/tarea0/K7/`.
