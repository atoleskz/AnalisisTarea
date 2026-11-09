# =============================================================================
#  Problema 5 (Integración numérica 21.8)
#  ∫_0^3 (5 + 3 cos x) dx  con fórmulas de Newton-Cotes cerradas y abiertas
#
#  Archivo: P05_newton_cotes.jl        (Tarea de Análisis Numérico 2027-1)
#  Uso:     julia P05_newton_cotes.jl
#  Salida:  tabla en consola, figuras/P05_*.png
# =============================================================================
using Printf
using Plots

f(x) = 5 + 3cos(x)
a, b = 0.0, 3.0
exacto = 5(b - a) + 3(sin(b) - sin(a))       # = 15 + 3 sen 3

# ---------------------- Fórmulas cerradas (una aplicación) -------------------
trapecio(f, a, b)   = (b - a) * (f(a) + f(b)) / 2
simpson13(f, a, b)  = (b - a) * (f(a) + 4f((a + b) / 2) + f(b)) / 6
function simpson38(f, a, b)
    h = (b - a) / 3
    return (b - a) * (f(a) + 3f(a + h) + 3f(a + 2h) + f(b)) / 8
end
function boole(f, a, b)
    h = (b - a) / 4
    x = [a + i * h for i in 0:4]
    return (b - a) * (7f(x[1]) + 32f(x[2]) + 12f(x[3]) + 32f(x[4]) + 7f(x[5])) / 90
end

# Aplicación múltiple de Simpson con n = 5 segmentos (n impar):
# Simpson 1/3 en los 2 primeros segmentos + Simpson 3/8 en los 3 últimos.
function simpson_multiple_n5(f, a, b)
    h = (b - a) / 5
    x = [a + i * h for i in 0:5]
    return simpson13(f, x[1], x[3]) + simpson38(f, x[3], x[6])
end

# ---------------------- Fórmulas abiertas (Newton-Cotes) ---------------------
# h = (b - a)/n_segmentos; los extremos a y b NO se evalúan.
function punto_medio(f, a, b)                 # 2 segmentos, 1 punto
    return (b - a) * f((a + b) / 2)
end
function abierta_3seg_2pts(f, a, b)           # (b - a)(f1 + f2)/2
    h = (b - a) / 3
    return (b - a) * (f(a + h) + f(a + 2h)) / 2
end
function abierta_4seg_3pts(f, a, b)           # (b - a)(2f1 - f2 + 2f3)/3
    h = (b - a) / 4
    return (b - a) * (2f(a + h) - f(a + 2h) + 2f(a + 3h)) / 3
end

# Términos de error teóricos (con f^(k)(ξ) acotada por max|f^(k)| = 3 en [0,3])
#  trapecio:      -(b-a)^3/12 f''          Simpson 1/3: -(b-a)^5/2880 f''''
#  Simpson 3/8:   -(b-a)^5/6480 f''''      Boole:       -2(b-a)^7/945 f^(6) / 4^7
#  punto medio:    (b-a)^3/24 f''          abierta 3-2:  3(b-a)^3/108 f''
#  abierta 4-3:    14(b-a)^5/45/4^5 f''''
L = b - a
cotas = Dict(
    "a) Trapecio (1 aplicación)"         => L^3 / 12 * 3,
    "b) Simpson 1/3"                     => L^5 / 2880 * 3,
    "c) Simpson 3/8"                     => L^5 / 6480 * 3,
    "d) Simpson múltiple, n = 5"         => (2(L/5))^5 / 2880 * 3 + (3(L/5))^5 / 6480 * 3,
    "e) Boole"                           => 8 / 945 * (L / 4)^7 * 3,
    "f) Punto medio"                     => L^3 / 24 * 3,
    "g) Abierta, 3 segmentos y 2 puntos" => 3 / 4 * (L / 3)^3 * 3,
    "h) Abierta, 4 segmentos y 3 puntos" => 14 / 45 * (L / 4)^5 * 3,
)

metodos = [
    ("a) Trapecio (1 aplicación)",         trapecio),
    ("b) Simpson 1/3",                     simpson13),
    ("c) Simpson 3/8",                     simpson38),
    ("d) Simpson múltiple, n = 5",         simpson_multiple_n5),
    ("e) Boole",                           boole),
    ("f) Punto medio",                     punto_medio),
    ("g) Abierta, 3 segmentos y 2 puntos", abierta_3seg_2pts),
    ("h) Abierta, 4 segmentos y 3 puntos", abierta_4seg_3pts),
]

@printf("Valor analítico: 15 + 3 sen(3) = %.10f\n\n", exacto)
@printf("%-36s %14s %14s %12s %12s\n", "método", "aproximación", "E_t", "ε_t (%)", "cota |E|")
resultados = Tuple{String,Float64,Float64}[]
for (nom, regla) in metodos
    I = regla(f, a, b)
    Et = exacto - I
    εt = abs(Et / exacto) * 100
    push!(resultados, (nom, I, εt))
    @printf("%-36s %14.8f %14.4e %12.6f %12.4e\n", nom, I, Et, εt, cotas[nom])
end

# Extra: convergencia de las reglas compuestas (trapecio y Simpson 1/3)
function trapecio_comp(f, a, b, n)
    h = (b - a) / n
    return h * (f(a) / 2 + sum(f(a + i * h) for i in 1:n-1; init = 0.0) + f(b) / 2)
end
function simpson_comp(f, a, b, n)
    h = (b - a) / n
    return h / 3 * (f(a) + f(b) + 4sum(f(a + i * h) for i in 1:2:n-1) + 2sum(f(a + i * h) for i in 2:2:n-2; init = 0.0))
end
ns = [2, 4, 8, 16, 32, 64, 128]
eT = [abs(exacto - trapecio_comp(f, a, b, n)) for n in ns]
eS = [abs(exacto - simpson_comp(f, a, b, n)) for n in ns]
println("\nExtra: reglas compuestas")
@printf("%5s %14s %14s\n", "n", "error trapecio", "error Simpson")
for i in eachindex(ns)
    @printf("%5d %14.4e %14.4e\n", ns[i], eT[i], eS[i])
end
@printf("Pendientes log-log: trapecio %.2f, Simpson %.2f\n",
        log(eT[end-1] / eT[1]) / log(ns[1] / ns[end-1]), log(eS[4] / eS[1]) / log(ns[1] / ns[4]))

# ------------------------------- Gráficas ------------------------------------
mkpath(joinpath(@__DIR__, "..", "figuras"))
etiquetas = [split(r[1], ")")[1] * ")" for r in resultados]
p1 = bar(etiquetas, [r[3] for r in resultados], yscale = :log10, legend = false,
         ylabel = "ε_t (%)", title = "21.8 Error relativo porcentual por método",
         color = [:firebrick, :steelblue, :steelblue, :steelblue, :seagreen, :darkorange, :darkorange, :darkorange])
xx = range(a, b, length = 300)
p2 = plot(xx, f.(xx), lw = 3, color = :black, label = "f(x) = 5 + 3cos x", xlabel = "x", title = "Integrando y nodos")
plot!(p2, [a, b], [f(a), f(b)], label = "trapecio", color = :firebrick, fill = (0, 0.12, :firebrick))
h4 = (b - a) / 4
scatter!(p2, [a + i * h4 for i in 0:4], f.([a + i * h4 for i in 0:4]), label = "nodos de Boole", color = :seagreen, ms = 6)
scatter!(p2, [a + i * h4 for i in 1:3], f.([a + i * h4 for i in 1:3]), label = "nodos abierta 4-3", color = :darkorange, marker = :x, ms = 8)
p3 = plot(ns, eT, xscale = :log2, yscale = :log10, marker = :circle, label = "trapecio compuesto",
          xlabel = "n", ylabel = "|E_t|", title = "Convergencia (extra)")
plot!(p3, ns, max.(eS, 1e-16), marker = :square, label = "Simpson 1/3 compuesto")
fig = plot(p1, p2, p3, layout = (1, 3), size = (1450, 430), margin = 6Plots.mm)
savefig(fig, joinpath(@__DIR__, "..", "figuras", "P05_newton_cotes.png"))
println("Figura guardada en figuras/P05_newton_cotes.png")
