# Tarea A.2 - diagnostico de OpenMP, ronda 1 (como pedir los nucleos)

Corridas cortas del caso base (NELM = 4), enviadas a la vez. Topologia: regla de Granado (c y c+36 = mismo nucleo).

| Variante | Peticion a Slurm | p | CPUs cobradas | Mascara (nucleos fisicos / CPUs logicos) | Iteraciones | t/iter (s) | CPU/real | Veredicto |
|---|---|---|---|---|---|---|---|---|
| A | `--ntasks=1 --cpus-per-task=p --hint=nomultithread` | 4 | 4 | 2 / 4 | - | - | - | NO PASA: mascara con 2 nucleos fisicos de 4; Slurm la corto (TIMEOUT): no termino las 4 iteraciones en 20 min |
| B | `--ntasks=1 --cpus-per-task=p --threads-per-core=1` | 4 | 4 | 2 / 4 | - | - | - | NO PASA: mascara con 2 nucleos fisicos de 4; Slurm la corto (TIMEOUT): no termino las 4 iteraciones en 20 min |
| C | `--ntasks=p --cpus-per-task=1 --hint=nomultithread + srun -n 1 -c p` | 4 | 4 | 2 / 4 | - | - | - | NO PASA: mascara con 2 nucleos fisicos de 4; Slurm la corto (TIMEOUT): no termino las 4 iteraciones en 20 min |
| D | `--ntasks=1 --cpus-per-task=2p --hint=multithread` | 4 | 8 | 4 / 8 | 4 | 130.2 | 4.89 | TERMINO: mascara de 4 nucleos fisicos, 4 iteraciones completas |
| B | `--ntasks=1 --cpus-per-task=p --threads-per-core=1` | 24 | 24 | 12 / 24 | - | - | - | NO PASA: mascara con 12 nucleos fisicos de 24; Slurm la corto (TIMEOUT): no termino las 4 iteraciones en 20 min |
| C | `--ntasks=p --cpus-per-task=1 --hint=nomultithread + srun -n 1 -c p` | 24 | 24 | 12 / 24 | - | - | - | NO PASA: mascara con 12 nucleos fisicos de 24; Slurm la corto (TIMEOUT): no termino las 4 iteraciones en 20 min |

Referencia: 1 proceso con 1 hilo (A.1, p = 1), primeras 4 iteraciones: 95.8 s por iteracion.

Las corridas que Slurm corto por tiempo no dejan OSZICAR ni OUTCAR: VASP compilado con gfortran los vacia al terminar. Su mascara se leyo de `vasp.err` (lineas de `srun --cpu-bind=verbose`).
