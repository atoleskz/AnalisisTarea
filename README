# Tarea de Análisis Numérico 2027-1: códigos en Julia

Seis problemas resueltos en Julia, sin bibliotecas de integración ni de interpolación.

| Archivo | Problema | Métodos |
|---|---|---|
| `codigos/P01_horner_error.jl` | Errores, ej. 11 | Horner, cota de error a posteriori, BigFloat |
| `codigos/P02_euler_rk4_adams.jl` | EDO 9.1, 9.22, 9.32, 9.33 | Euler, RK4, predictor-corrector de Adams (orden 3 y 4) |
| `codigos/P03_proyectil_rk.jl` | EDO 9.18 | RK2 (traducción del FORTRAN), RK3, regla del trapecio |
| `codigos/P04_disparo_rk4.jl` | EDO 9.31 | Disparo, RK4 para sistemas, bisección + secante |
| `codigos/P05_newton_cotes.jl` | Integración 21.8 | Newton-Cotes cerradas y abiertas, Boole |
| `codigos/P06_splines_acciones.jl` | Splines, computadora 14 | Diferencias divididas de Newton, spline cúbica natural (Thomas) |

## Cómo correrlos

```bash
julia -e 'using Pkg; Pkg.add("Plots")'   # una sola vez
julia codigos/P01_horner_error.jl
julia codigos/P02_euler_rk4_adams.jl
# ...
```

Cada programa imprime sus tablas en la consola y guarda sus gráficas en `figuras/`.
`P06` lee `datos/aapl_101_cierres.csv`, que tiene 101 cierres diarios de AAPL del 17-feb-2015
al 10-jul-2015 (del conjunto `plotly/datasets`, originalmente de Yahoo Finance).

En Colab (entorno de Julia) se pueden pegar los archivos en celdas. Las rutas usan `@__DIR__`;
si se corre desde una celda, hay que crear las carpetas `figuras/` y `datos/` junto al notebook.

## Nota sobre el problema 9.31
Las constantes `A, k, P, Q` están en el bloque `CONSTANTES` de `P04_disparo_rk4.jl`. Deben ser
las del ejemplo 9.14 del libro de Nakamura; si son otras, sólo hay que cambiarlas ahí.
