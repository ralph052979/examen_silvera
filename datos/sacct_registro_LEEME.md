# Registro de consumo: qué jobs son del examen

`sacct_registro_examen.txt` es la salida del comando que pide el enunciado, corrido en Granado:

```
sacct -u $USER -S 2026-09-07 -X --format=JobID,JobName,Partition,AllocCPUS,Elapsed,CPUTimeRAW,MaxRSS,State
```

La cuenta `rsilvera` es también la de mi investigación. Como el comando empieza el 7 de septiembre, el
registro contiene dos grupos de jobs:

| Jobs | JobID | Cantidad | h-núcleo (CPUTimeRAW / 3600) |
|---|---|---|---|
| **Del examen** (17 al 22 de noviembre) | ≥ 528789 | 166 | **148.5 de 250** |
| Investigación anterior al examen | < 528789 | 370 | 11 682.7 (fuera del presupuesto del examen) |

* El primer job del examen es el 528789 (`calib_8p`, enviado el 2026-11-17 a las 19:41). Los 166 jobs desde
  ahí son todos del examen, incluidas las pruebas y los intentos previos del 17 y 18 de noviembre. El consumo
  por parte está en `parteD_tabla.md` y `parteD_consumo.csv`.
* Los 370 jobs anteriores son cálculos de investigación sobre WSe₂: relajaciones, SCF y NSCF, bandas,
  convergencia en ENCUT y en k-points, fonones y pruebas (`WSe2_relax`, `RSA-scf`, `RSA-ncf`, `WSe2_bands`,
  `conv_E_*`, `conv_K_*`, `fon_*`, ...).

La separación se comprueba desde el propio registro (es lo que hace `slurm/notas_finales.sh`):

```
awk '$1~/^[0-9]/ {h=$6/3600; if ($1+0<528789) {n1++; h1+=h} else {n2++; h2+=h}} END {printf "antes de 528789: %d jobs, %.1f h\ndesde 528789: %d jobs, %.1f h\n", n1,h1,n2,h2}' datos/sacct_registro_examen.txt
```
