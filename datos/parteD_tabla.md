# Parte D - consumo, colas y costo por calculo

Consumo en phys_fast y phys_general desde 2026-11-16 (sacct, CPUTimeRAW): **148.53 h-nucleo** de 250.

| parte | jobs | h-nucleo |
|---|---|---|
| Tarea 0: calibracion final (K1-K7) | 7 | 2.00 |
| A.1: escalamiento fuerte MPI | 21 | 21.59 |
| A.2: diagnostico de OpenMP (2 rondas) | 18 | 33.95 |
| A.2: barrido OpenMP (4 iteraciones) | 21 | 30.94 |
| B.1: escalamiento debil | 21 | 15.55 |
| B.2: hibrido MPI x OpenMP | 19 | 23.04 |
| C.2: serie de ENCUT | 3 | 2.37 |
| Pruebas e intentos previos (17-18 nov) | 56 | 19.09 |
| **total** | 166 | **148.53** |

Espera en cola de los 143 jobs del examen (Start - Eligible): mediana 0 s, percentil 90 0 s, maximo 12.7 min; 98 % empezaron en menos de 1 minuto.

| configuracion | nucleos | t (s) | h-nucleo por calculo |
|---|---|---|---|
| 1 procesos MPI, KPAR = 1 (A.1) | 1 | 1871.6 +/- 1.7 | 0.52 |
| 2 procesos MPI, KPAR = 1 (A.1) | 2 | 1017.7 +/- 1.1 | 0.57 |
| 4 procesos MPI, KPAR = 1 (A.1) | 4 | 579.8 +/- 4.0 | 0.64 |
| 8 procesos MPI, KPAR = 1 (A.1) | 8 | 378.1 +/- 2.2 | 0.84 |
| 12 procesos MPI, KPAR = 1 (A.1) | 12 | 329.4 +/- 4.1 | 1.10 |
| 16 procesos MPI, KPAR = 1 (A.1) | 16 | 328.8 +/- 3.8 | 1.46 |
| 24 procesos MPI, KPAR = 1 (A.1) | 24 | 288.5 +/- 6.0 | 1.92 |
| 12 procesos MPI, KPAR = 12 (B.1, malla del caso base) | 12 | 235.8 +/- 0.3 | 0.79 |
