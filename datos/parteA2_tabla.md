# Tarea A.2 - escalamiento fuerte con OpenMP

Caso base: WSe2 bulk 2x2x1 con V_Se (23 atomos), ENCUT = 400 eV, Gamma 10x10x3 (36 k-points), NBANDS = 144. 1 proceso con p hilos OpenMP (OMP_PLACES=cores, OMP_PROC_BIND=close), un nucleo fisico por hilo; nucleos pedidos como en la A.1 y paso de srun elegido por su mascara real; OPENBLAS_NUM_THREADS=1; carga fija de 4 iteraciones SCF (NELM = 4). t = tiempo de pared de `srun` (arranque + 4 iteraciones + fuerzas y salida), media +/- desviacion estandar de 3 repeticiones; t(1) con la misma carga.

| p | n | t (s) | S(p) | E(p) | t/iteracion (s) | S por iteracion | E por iteracion | memoria (MB) |
|---|---|---|---|---|---|---|---|---|
| 1 | 3 | 590.8 +/- 5.1 | 1.00 +/- 0.01 | 1.000 +/- 0.012 | 95.5 +/- 0.7 | 1.00 +/- 0.01 | 1.000 | 933 |
| 2 | 3 | 476.1 +/- 4.4 | 1.24 +/- 0.02 | 0.620 +/- 0.008 | 76.8 +/- 1.1 | 1.24 +/- 0.02 | 0.622 | 936 |
| 4 | 3 | 426.4 +/- 7.2 | 1.39 +/- 0.03 | 0.346 +/- 0.007 | 70.0 +/- 0.6 | 1.36 +/- 0.01 | 0.341 | 943 |
| 8 | 3 | 433.9 +/- 5.5 | 1.36 +/- 0.02 | 0.170 +/- 0.003 | 75.7 +/- 0.8 | 1.26 +/- 0.02 | 0.158 | 965 |
| 12 | 3 | 503.7 +/- 5.2 | 1.17 +/- 0.02 | 0.098 +/- 0.001 | 94.3 +/- 1.2 | 1.01 +/- 0.01 | 0.084 | 971 |
| 16 | 3 | 516.6 +/- 0.8 | 1.14 +/- 0.01 | 0.071 +/- 0.001 | 97.0 +/- 0.2 | 0.98 +/- 0.01 | 0.062 | 1010 |
| 24 | 3 | 644.4 +/- 2.6 | 0.92 +/- 0.01 | 0.038 +/- 0.000 | 129.4 +/- 0.7 | 0.74 +/- 0.01 | 0.031 | 1023 |

Por iteracion SCF, MPI (A.1: p procesos con 1 hilo, primeras 4 iteraciones) frente a OpenMP (A.2: 1 proceso con p hilos):

| p | t/iter MPI (s) | S_iter MPI | E_iter MPI | t/iter OpenMP (s) | S_iter OpenMP | E_iter OpenMP |
|---|---|---|---|---|---|---|
| 1 | 95.8 | 1.00 +/- 0.00 | 1.000 | 95.5 | 1.00 +/- 0.01 | 1.000 |
| 2 | 52.4 | 1.83 +/- 0.01 | 0.913 | 76.8 | 1.24 +/- 0.02 | 0.622 |
| 4 | 29.6 | 3.23 +/- 0.03 | 0.808 | 70.0 | 1.36 +/- 0.01 | 0.341 |
| 8 | 19.4 | 4.95 +/- 0.05 | 0.618 | 75.7 | 1.26 +/- 0.02 | 0.158 |
| 12 | 16.8 | 5.71 +/- 0.08 | 0.476 | 94.3 | 1.01 +/- 0.01 | 0.084 |
| 16 | 16.9 | 5.65 +/- 0.09 | 0.353 | 97.0 | 0.98 +/- 0.01 | 0.062 |
| 24 | 14.3 | 6.71 +/- 0.13 | 0.280 | 129.4 | 0.74 +/- 0.01 | 0.031 |
