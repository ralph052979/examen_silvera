# Parte I — Paralelismo: modelos y métricas

Rafael Silvera Altamiranda. Los números se pueden verificar con
`python3 teoria/parte1_calculos.py`, que imprime todos los valores de esta página.

---

## Problema 1. Leyes de Amdahl y Gustafson

### (a) Deducción de la ley de Amdahl

Sea $t(1)$ el tiempo de un programa en un procesador. Una fracción $f$ de ese tiempo es serial
(no se puede repartir) y la fracción $1-f$ es paralelizable. Con $p$ procesadores, en el mejor
caso la parte paralelizable se reparte en partes iguales y sin costo adicional (sin comunicación
ni desbalance), mientras que la parte serial no cambia:

$$
t(p) = f\,t(1) + \frac{(1-f)\,t(1)}{p}.
$$

La aceleración es el cociente de los tiempos, y $t(1)$ se cancela:

$$
S(p) = \frac{t(1)}{t(p)} = \frac{1}{f + \dfrac{1-f}{p}}.
$$

Como $(1-f)/p \ge 0$, se cumple $S(p) < 1/f$ para todo $p$, y $S(p) \to 1/f$ cuando
$p \to \infty$. Ese es el límite de Amdahl: aunque haya infinitos procesadores, el tiempo nunca
baja de $f\,t(1)$. Cualquier sobrecosto real (comunicación, sincronización, desbalance) solo
puede bajar $S(p)$ por debajo de esta cota.

### (b) $f = 0.05$

$$
S(24) = \frac{1}{0.05 + 0.95/24} = \frac{1}{0.0895833} = 11.16,
\qquad E(24) = \frac{S(24)}{24} = 0.465,
\qquad S(\infty) = \frac{1}{0.05} = 20.
$$

Con solo un 5 % de código serial, 24 núcleos rinden menos de la mitad de lo ideal, y ningún
número de núcleos pasa de 20.

### (c) Ley de Gustafson

Gustafson (1988) parte de la ejecución en paralelo: con $p$ procesadores el tiempo se normaliza a
$s + (1-s) = 1$, donde $s$ es la fracción serial **medida en la corrida paralela**. Si el problema
se agranda con $p$ (más datos, más k-points, una supercelda mayor), el mismo trabajo en un solo
procesador tomaría $s + (1-s)\,p$, porque la parte paralela de cada procesador se haría una tras
otra. La aceleración escalada es

$$
S_G(p) = \frac{s + (1-s)\,p}{s + (1-s)} = s + (1-s)\,p = p - s\,(p-1).
$$

Es lineal en $p$, con pendiente $1-s$; por ejemplo, con $s = 0.05$, $S_G(24) = 22.85$. Es más
optimista porque supone que lo que crece con $p$ es solo la parte paralela: el trabajo serial
(leer la entrada, inicializar, escribir resultados) queda fijo, así que su peso relativo **baja**
al agrandar el problema. En Amdahl ocurre lo contrario: el problema es fijo y, al repartir la parte
paralela entre más procesadores, el peso relativo de la parte serial **sube**. Las dos fracciones
no son la misma: $s$ se mide en la corrida paralela y $f$ en la serial; para el mismo problema,
$f = s/[\,s + (1-s)\,p\,]$.

### (d) Tipo de escalamiento

* **Amdahl → escalamiento fuerte:** tamaño total del problema fijo; se mide $S(p) = t(1)/t(p)$.
  Es la Parte II‑A del examen.
* **Gustafson → escalamiento débil:** trabajo por procesador fijo, el problema crece con $p$; se
  mide la eficiencia débil $E_w(p) = t(1)/t(p)$ con el problema escalado. Es la Tarea B.1.

---

## Problema 2. Memoria compartida (OpenMP) frente a memoria distribuida (MPI)

### (a) Qué se comparte y qué no

| | OpenMP (memoria compartida) | MPI (memoria distribuida) |
|---|---|---|
| Unidad de ejecución | hilos de **un mismo proceso** | **procesos** independientes |
| Se comparte | todo el espacio de direcciones: variables globales y estáticas, el *heap* (arreglos dinámicos), el código y los archivos abiertos | nada en memoria: cada proceso tiene su propio espacio de direcciones (salvo, de forma explícita, las ventanas de memoria compartida de MPI‑3) |
| Es privado | la pila de cada hilo (variables locales), sus registros y contador de programa, y las variables declaradas `private`/`threadprivate` | todo; lo que varios procesos necesitan se **replica** o se **reparte** explícitamente |
| Comunicación | implícita: leer y escribir variables compartidas | explícita: mensajes punto a punto (`MPI_Send`/`MPI_Recv`) y colectivas (`MPI_Allreduce`, …) |
| Riesgos y costos típicos | condiciones de carrera, sincronización (barreras, `critical`, `atomic`, `reduction`), coherencia de caché, *false sharing*, efectos NUMA | latencia y ancho de banda de los mensajes, desbalance de carga, memoria replicada en cada proceso |

En este examen se ve en la memoria: con 24 procesos MPI la suma de la memoria de VASP es 4 veces
la de 1 proceso con 24 hilos (3.95 frente a 1.00 GB en la B.2), porque cada proceso MPI guarda su
propia copia de las mallas y los proyectores, mientras que los hilos comparten una sola copia.

### (b) Límite natural de OpenMP en un nodo de Granado

Los hilos necesitan memoria compartida, que solo existe dentro de un nodo: un proceso OpenMP no puede
usar más que los **36 núcleos físicos** del nodo (72 hilos de hardware con *hyperthreading*). En la
práctica el límite útil es de 36 hilos, uno por núcleo físico, porque los dos hilos de hardware de
un núcleo comparten sus unidades de cálculo, cachés y ancho de banda. En la A.2 de este examen, cuando
Slurm puso 4 hilos en 2 núcleos físicos (2 hilos de hardware por núcleo), las 4 iteraciones no
terminaron en 20 min, más lento que con un solo hilo (95.8 s por iteración); con 4 núcleos físicos y
un hilo por núcleo, cada iteración tardó 69 s. Además, con más de 18 hilos el proceso ocupa los dos sockets y parte de
sus accesos a memoria van al banco NUMA del otro socket, así que el punto natural para un proceso con
hilos es un socket (hasta 18 hilos).

### (c) Por qué MPI escala a varios nodos y OpenMP no

MPI no necesita un espacio de direcciones común: cada proceso trabaja con su memoria local y todo
intercambio es un mensaje, que puede viajar por la red de interconexión entre nodos igual que por
memoria compartida dentro de un nodo. OpenMP supone que todos los hilos pueden leer y escribir las
mismas direcciones de memoria con coherencia de caché, y eso lo da el hardware solo dentro de un nodo;
la memoria de otro nodo no es direccionable por los hilos. Por eso el modelo habitual en un clúster
es híbrido: MPI entre nodos (o entre sockets) y, opcionalmente, OpenMP dentro de cada proceso.

---

## Problema 3. Eficiencia paralela con tiempos medidos

La eficiencia paralela es $E(p) = S(p)/p = t(1)/[\,p\,t(p)\,]$: la fracción de los $p$ núcleos que
se aprovecha (1 = ideal).

| $p$ | $t(p)$ (s) | $S(p) = t(1)/t(p)$ | $E(p) = S(p)/p$ | Karp–Flatt $e(p) = \dfrac{1/S - 1/p}{1 - 1/p}$ |
|---|---|---|---|---|
| 1 | 100 | 1.000 | 1.000 | — |
| 4 | 30 | 3.333 | 0.833 | 0.0667 |
| 16 | 11 | 9.091 | 0.568 | 0.0507 |
| 24 | 9 | 11.111 | **0.463** | 0.0504 |

**¿Desde cuántos núcleos cae la eficiencia por debajo del 50 %?** Entre los puntos medidos, desde
$p = 24$: $E(16) = 0.568$ y $E(24) = 0.463$. Con la $f$ estimada abajo, Amdahl da
$E(p) = 1/(f\,p + 1 - f) = 0.5$ en $p = (1 + f)/f \approx 19$–$20$ núcleos.

**Fracción serial efectiva.** Despejando $f$ de la ley de Amdahl para cada punto se obtiene la
fracción de Karp–Flatt $e(p)$ de la tabla: 0.067, 0.051 y 0.050. Ajustando por mínimos cuadrados
$t(p) = a + b/p$ a los cuatro tiempos (con $a$ la parte serial y $b$ la paralelizable):

$$
a = 5.42\ \text{s},\qquad b = 94.77\ \text{s},\qquad f = \frac{a}{a + b} = 0.054.
$$

El ajuste lineal de $1/S - 1/p = f\,(1 - 1/p)$ da el mismo valor, $f = 0.054$. Es decir, $f \approx 5\,\%$,
el caso del Problema 1: Amdahl con $f = 0.05$ predice $S(24) = 11.16$ y se midió $11.11$. Que
$e(p)$ no crezca con $p$ indica que, en este ejemplo, el límite es la fracción serial y no un
sobrecosto que aumente con $p$ (como la comunicación). En el caso de VASP de la Parte II sí crece
(ver la Tarea A.3 del reporte).
