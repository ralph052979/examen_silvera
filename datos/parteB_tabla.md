# Parte B - escalamiento debil e hibrido

## B.1 Escalamiento debil

p procesos MPI con KPAR = p (un grupo de k-points por proceso) y una malla Gamma con 3p k-points irreducibles: cada proceso calcula 3 k-points con las 144 bandas. ENCUT = 400 eV, SCF completo (EDIFF = 1E-6). Ew(p) = t(1)/t(p). Memoria: VmHWM de cada proceso de VASP, leida cada 3 s.

| p | malla | NKPTS | n | t (s) | iter. SCF | t/iter (s) | Ew | Ew por iter. | mem/proceso (MB) | mem total (GB) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 3x3x1 | 3 | 3 | 164.5 +/- 0.5 | 18/18/18 | 8.03 +/- 0.03 | 1.000 +/- 0.005 | 1.000 +/- 0.005 | 183 | 0.18 |
| 2 | 3x3x2 | 6 | 3 | 169.4 +/- 0.5 | 18/18/18 | 8.27 +/- 0.03 | 0.971 +/- 0.004 | 0.971 +/- 0.005 | 250 | 0.49 |
| 4 | 5x5x3 | 12 | 3 | 180.0 +/- 0.8 | 18/18/18 | 8.74 +/- 0.04 | 0.914 +/- 0.005 | 0.918 +/- 0.006 | 388 | 1.51 |
| 8 | 6x6x4 | 24 | 3 | 205.8 +/- 0.3 | 18/18/18 | 9.94 +/- 0.01 | 0.799 +/- 0.003 | 0.807 +/- 0.003 | 662 | 5.16 |
| 12 | 10x10x3 | 36 | 3 | 235.8 +/- 0.3 | 18/18/18 | 11.34 +/- 0.02 | 0.698 +/- 0.003 | 0.708 +/- 0.003 | 936 | 10.95 |
| 16 | 14x14x2 | 48 | 3 | 284.4 +/- 0.7 | 18/18/18 | 13.62 +/- 0.04 | 0.578 +/- 0.002 | 0.589 +/- 0.003 | 1211 | 18.88 |
| 24 | 10x10x6 | 72 | 3 | 327.3 +/- 1.0 | 18/18/18 | 15.42 +/- 0.04 | 0.503 +/- 0.002 | 0.520 +/- 0.002 | 1765 | 41.33 |

## B.2 Hibrido MPI x OpenMP con 24 nucleos

Carga fija de 4 iteraciones SCF del caso base (NELM = 4); t(1) = 590.8 +/- 5.1 s es la p = 1 de la A.2 con la misma carga. (1,24) incluye las corridas p = 24 de la A.2. Cada proceso con T nucleos fisicos propios (mascara verificada), OMP_PLACES=cores, OMP_PROC_BIND=close, OPENBLAS_NUM_THREADS=1.

| N_MPI x N_OMP | n | t (s) | t/iter (s) | S = t(1)/t | E = S/24 | mem/proceso (MB) | mem total (GB) | procesos en 2 sockets |
|---|---|---|---|---|---|---|---|---|
| 24 x 1 | 3 (3 B.2) | 89.3 +/- 1.3 | 14.30 +/- 0.31 | 6.62 +/- 0.11 | 0.276 | 174 | 3.95 | 0/0/0 |
| 12 x 2 | 3 (3 B.2) | 92.0 +/- 0.8 | 14.59 +/- 0.21 | 6.42 +/- 0.08 | 0.267 | 206 | 2.36 | 0/0/0 |
| 8 x 3 | 3 (3 B.2) | 107.1 +/- 0.4 | 17.14 +/- 0.04 | 5.51 +/- 0.05 | 0.230 | 239 | 1.84 | 0/0/0 |
| 6 x 4 | 3 (3 B.2) | 123.9 +/- 0.9 | 20.64 +/- 0.23 | 4.77 +/- 0.05 | 0.199 | 279 | 1.62 | 1/1/1 |
| 4 x 6 | 3 (3 B.2) | 152.5 +/- 1.6 | 25.83 +/- 0.38 | 3.87 +/- 0.05 | 0.161 | 348 | 1.35 | 0/0/0 |
| 2 x 12 | 3 (3 B.2) | 311.0 +/- 1.7 | 58.39 +/- 0.25 | 1.90 +/- 0.02 | 0.079 | 572 | 1.11 | 1/1/1 |
| 1 x 24 | 4 (1 B.2+3 A.2) | 644.6 +/- 2.2 | 129.49 +/- 0.56 | 0.92 +/- 0.01 | 0.038 | 1020 | 1.00 | 1 |

Combinacion optima (menor t medio): **24 x 1**, t = 89.3 +/- 1.3 s.
