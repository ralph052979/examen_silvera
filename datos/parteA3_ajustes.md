# Tarea A.3 - ajustes de Amdahl

Amdahl: 1/S = f + (1-f)/p. Con sobrecosto: 1/S = f + (1-f)/p + c(p-1). Minimos cuadrados ponderados sobre 1/S - 1/p (p > 1); incertidumbres de 1 sigma, agrandadas por sqrt(chi2 reducido) si este pasa de 1.

| serie | f (Amdahl) | chi2 red. | S(24) Amdahl | S max = 1/f | f (con sobrecosto) | c | chi2 red. | p* | S(p*) |
|---|---|---|---|---|---|---|---|---|---|
| MPI corrida | 0.0948 +/- 0.0059 | 53.36 | 7.55 | 10.5 | 0.0787 +/- 0.0057 | 1.89e-03 +/- 5.5e-04 | 16.72 | 22.1 | 6.24 |
| OpenMP corrida | 0.8522 +/- 0.0673 | 189.22 | 1.16 | 1.2 | 0.5329 +/- 0.0308 | 2.25e-02 +/- 2.0e-03 | 7.01 | 4.6 | 1.40 |
| MPI por iteracion | 0.1010 +/- 0.0063 | 24.32 | 7.22 | 9.9 | 0.0796 +/- 0.0094 | 1.70e-03 +/- 6.6e-04 | 11.51 | 23.3 | 6.37 |
| OpenMP por iteracion | 0.9351 +/- 0.1062 | 507.80 | 1.07 | 1.1 | 0.4978 +/- 0.0377 | 3.47e-02 +/- 2.6e-03 | 14.32 | 3.8 | 1.38 |

Fraccion serial de Karp-Flatt e(p) = (1/S - 1/p)/(1 - 1/p):

| serie | p = 2 | p = 4 | p = 8 | p = 12 | p = 16 | p = 24 |
|---|---|---|---|---|---|---|
| MPI corrida | 0.088 | 0.080 | 0.088 | 0.101 | 0.121 | 0.117 |
| OpenMP corrida | 0.612 | 0.629 | 0.696 | 0.839 | 0.866 | 1.095 |
| MPI por iteracion | 0.095 | 0.079 | 0.088 | 0.100 | 0.122 | 0.112 |
| OpenMP por iteracion | 0.608 | 0.644 | 0.763 | 0.987 | 1.017 | 1.371 |
