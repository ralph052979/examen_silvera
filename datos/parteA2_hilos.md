# Tarea A.2 - diagnostico de OpenMP, ronda 2 (hilos de OpenBLAS y espera de OpenMP)

Intento 1, variante B (`--ntasks=1 --cpus-per-task=p --threads-per-core=1`): E1 p=4 mascara 2/4, afinidad; E3 p=4 mascara 2/4, afinidad; E2 p=4 mascara 2/4, afinidad; E0 p=4 mascara 2/4, afinidad; E1 p=24 omitida; E3 p=24 omitida. La mascara se reviso al arrancar y las corridas se detuvieron en segundos.

Intento 2: 1 proceso con p hilos; nucleos pedidos como en la A.1 (`--ntasks=p --cpus-per-task=1 --hint=nomultithread`); antes de VASP se prueban varias formas de lanzar el paso con srun y se usa la primera cuya mascara real tiene p nucleos fisicos (`geometria_pasos.txt` de cada corrida). NELM = 4, en una sola cadena.

| Entorno | Variables | p | Paso de srun | Mascara (nucleos) | Hilos activos (mediana) | Uso de CPU | Iteraciones | t/iter (s) | S por iteracion | Veredicto |
|---|---|---|---|---|---|---|---|---|---|---|
| E1 | OPENBLAS_NUM_THREADS=1 | 4 | hint | 4 | 4 | 3.92 | 4 | 69.8 | 1.37 | PASA: 4 hilos activos para 4 nucleos |
| E3 | OPENBLAS_NUM_THREADS=1 + OMP_WAIT_POLICY=PASSIVE | 4 | hint | 4 | 1 | 1.85 | 4 | 103.2 | 0.93 | LENTA: mas lenta que con 1 hilo (1 hilos activos para 4 nucleos) |
| E2 | OMP_WAIT_POLICY=PASSIVE | 4 | hint | 4 | 1 | 1.85 | 4 | 103.9 | 0.92 | LENTA: mas lenta que con 1 hilo (1 hilos activos para 4 nucleos) |
| E0 | de la ronda 1 (control) | 4 | hint | 4 | 4 | 3.91 | 4 | 69.1 | 1.39 | PASA: 4 hilos activos para 4 nucleos |
| E1 | OPENBLAS_NUM_THREADS=1 | 24 | hint | 24 | 24 | 23.47 | 4 | 129.6 | 0.74 | LENTA: mas lenta que con 1 hilo (24 hilos activos para 24 nucleos) |
| E3 | OPENBLAS_NUM_THREADS=1 + OMP_WAIT_POLICY=PASSIVE | 24 | hint | 24 | 1 | 2.93 | 3 | 216.2 | 0.44 | NO PASA: limite de tiempo (3 iteraciones SCF) |

Referencia con 1 hilo (A.1, p = 1, primeras 4 iteraciones): 95.8 s por iteracion.
OpenBLAS: OpenBLAS 0.3.17 NO_AFFINITY SKYL SINGLE_THREADEDOPENBLAS_VERBOSE; con OpenMP (enlaza libgomp): no.
Pila de los hilos (gdb) en E0 p=4: gdb no esta instalado en el nodo.

**Conclusion.** E1 (`OPENBLAS_NUM_THREADS=1`, espera de libgomp por defecto) es correcto con p = 4 y 24 (mascara de p nucleos fisicos, p hilos activos), pero OpenMP casi no acelera este caso: S por iteracion = 1.37 con 4 hilos y 0.74 con 24. La espera pasiva es mas lenta. El barrido de la A.2 usa E1 con una carga fija de 4 iteraciones SCF (README, seccion 6).
