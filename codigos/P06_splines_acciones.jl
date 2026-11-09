# =============================================================================
#  Problema 6 (Splines cúbicos, problema de computadora 14)
#  101 precios diarios de cierre de una acción: polinomio interpolante de
#  grado 20 vs. spline cúbica natural con nodos 0:5:100.
#
#  Datos: datos/aapl_101_cierres.csv  (Apple Inc., NASDAQ: AAPL, cierres del
#         17-feb-2015 al 10-jul-2015, conjunto público "finance-charts-apple"
#         de plotly/datasets, originalmente de Yahoo Finance).
#         Para usar otra acción basta reemplazar el CSV (columnas: dia,fecha,cierre).
#
#  Archivo: P06_splines_acciones.jl    (Tarea de Análisis Numérico 2027-1)
#  Uso:     julia P06_splines_acciones.jl
#  Salida:  tablas en consola, figuras/P06_*.png
# =============================================================================
using Printf
using Plots

# ------------------------------ Lectura de datos -----------------------------
ruta = joinpath(@__DIR__, "..", "datos", "aapl_101_cierres.csv")
lineas = readlines(ruta)[2:end]
dias    = [parse(Float64, split(l, ",")[1]) for l in lineas]
fechas  = [split(l, ",")[2] for l in lineas]
precios = [parse(Float64, split(l, ",")[3]) for l in lineas]
@assert length(precios) == 101 "Se requieren 101 precios"
println("Acción: AAPL, del $(fechas[1]) al $(fechas[end]) ($(length(precios)) cierres)")

x0 = collect(0.0:5.0:100.0)            # nodos de interpolación (21 puntos)
y0 = precios[1:5:101]                  # precios en los días 0, 5, ..., 100
x  = collect(0.0:1.0:100.0)            # puntos de evaluación

# ------------------- (a) Polinomio de interpolación de grado 20 --------------
# Diferencias divididas de Newton. Se hace en Float64 y en BigFloat (256 bits)
# para separar el error de interpolación del error de redondeo.
function difdiv(xn, yn)
    n = length(xn)
    c = copy(yn)
    for j in 2:n, i in n:-1:j
        c[i] = (c[i] - c[i-1]) / (xn[i] - xn[i-j+1])
    end
    return c
end
function newton_eval(c, xn, t)
    p = c[end]
    for i in length(c)-1:-1:1
        p = p * (t - xn[i]) + c[i]
    end
    return p
end

c64 = difdiv(x0, y0)
cbig = difdiv(BigFloat.(x0), BigFloat.(y0))
P64 = [newton_eval(c64, x0, t) for t in x]
Pbig = Float64.([newton_eval(cbig, BigFloat.(x0), big(t)) for t in x])

err_P = abs.(Pbig .- precios)
iP = argmax(err_P)
@printf("\n(a) Polinomio de grado 20\n")
@printf("    máximo error |P20(x) - precio| = %.4f USD  en el día %d (%s)\n", err_P[iP], Int(x[iP]), fechas[iP])
@printf("    máximo error en [10, 90]       = %.4f USD\n", maximum(err_P[11:91]))
@printf("    diferencia Float64 vs BigFloat  = %.3e USD (redondeo)\n", maximum(abs.(P64 .- Pbig)))
@printf("    rango de P20 en [0,100]: [%.2f, %.2f]  vs. precios en [%.2f, %.2f]\n",
        minimum(Pbig), maximum(Pbig), minimum(precios), maximum(precios))
@printf("    ¿pasa por los nodos? max |P20(x_i) - y_i| = %.3e\n", maximum(abs.(Pbig[1:5:101] .- y0)))

# ------------------- (b) Spline cúbica natural --------------------------------
# S_i(x) = a_i + b_i (x - x_i) + c_i (x - x_i)^2 + d_i (x - x_i)^3,  c_1 = c_n = 0
# Sistema tridiagonal para c, resuelto con el algoritmo de Thomas.
function spline_natural(xn, yn)
    n = length(xn)
    h = diff(xn)
    # sistema para c_2 ... c_{n-1}
    m = n - 2
    sub = [h[i] for i in 2:m]                      # subdiagonal
    dia = [2(h[i] + h[i+1]) for i in 1:m]          # diagonal
    sup = [h[i+1] for i in 1:m-1]                  # superdiagonal
    rhs = [3((yn[i+2] - yn[i+1]) / h[i+1] - (yn[i+1] - yn[i]) / h[i]) for i in 1:m]
    # Thomas
    for i in 2:m
        w = sub[i-1] / dia[i-1]
        dia[i] -= w * sup[i-1]
        rhs[i] -= w * rhs[i-1]
    end
    cint = zeros(m)
    cint[m] = rhs[m] / dia[m]
    for i in m-1:-1:1
        cint[i] = (rhs[i] - sup[i] * cint[i+1]) / dia[i]
    end
    c = vcat(0.0, cint, 0.0)
    a = yn[1:end-1]
    b = [(yn[i+1] - yn[i]) / h[i] - h[i] * (2c[i] + c[i+1]) / 3 for i in 1:n-1]
    d = [(c[i+1] - c[i]) / (3h[i]) for i in 1:n-1]
    return a, b, c[1:end-1], d
end
function spline_eval(coef, xn, t)
    a, b, c, d = coef
    i = clamp(searchsortedlast(xn, t), 1, length(xn) - 1)
    s = t - xn[i]
    return a[i] + b[i] * s + c[i] * s^2 + d[i] * s^3
end

coef = spline_natural(x0, y0)
S = [spline_eval(coef, x0, t) for t in x]
err_S = abs.(S .- precios)
iS = argmax(err_S)
@printf("\n(b) Spline cúbica natural (nodos 0:5:100)\n")
@printf("    máximo error |S(x) - precio| = %.4f USD  en el día %d (%s)\n", err_S[iS], Int(x[iS]), fechas[iS])
@printf("    máximo error en [10, 90]     = %.4f USD\n", maximum(err_S[11:91]))
@printf("    rango de S en [0,100]: [%.2f, %.2f]\n", minimum(S), maximum(S))
# verificación de continuidad de S, S', S'' en los nodos interiores
a_, b_, c_, d_ = coef
h5 = 5.0
saltoS1 = maximum(abs(b_[i] + 2c_[i] * h5 + 3d_[i] * h5^2 - b_[i+1]) for i in 1:19)
saltoS2 = maximum(abs(2c_[i] + 6d_[i] * h5 - 2c_[i+1]) for i in 1:19)
@printf("    continuidad en nodos: max salto S' = %.2e, max salto S'' = %.2e\n", saltoS1, saltoS2)

# ------------------- (c) Comparación ----------------------------------------
rms(v) = sqrt(sum(v .^ 2) / length(v))
println("\n(c) Comparación")
@printf("%-26s %12s %12s %12s\n", "", "max error", "RMS", "max [10,90]")
@printf("%-26s %12.4f %12.4f %12.4f\n", "Polinomio grado 20", maximum(err_P), rms(err_P), maximum(err_P[11:91]))
@printf("%-26s %12.4f %12.4f %12.4f\n", "Spline cúbica natural", maximum(err_S), rms(err_S), maximum(err_S[11:91]))

# ------------------------------- Gráficas ------------------------------------
mkpath(joinpath(@__DIR__, "..", "figuras"))
xf = range(0, 100, length = 2001)
Pf = Float64.([newton_eval(cbig, BigFloat.(x0), big(t)) for t in xf])
Sf = [spline_eval(coef, x0, t) for t in xf]

p1 = plot(x, precios, color = :gray, lw = 1.5, marker = :circle, ms = 2, label = "cierres diarios",
          xlabel = "día", ylabel = "USD", title = "(a) Polinomio de interpolación de grado 20")
plot!(p1, xf, Pf, color = :firebrick, lw = 2, label = "P20(x)")
scatter!(p1, x0, y0, color = :black, ms = 4, label = "nodos 0:5:100")
p1z = plot(p1, ylims = (minimum(precios) - 5, maximum(precios) + 5), title = "(a) detalle (escala de los datos)")
p2 = plot(x, precios, color = :gray, lw = 1.5, marker = :circle, ms = 2, label = "cierres diarios",
          xlabel = "día", ylabel = "USD", title = "(b) Spline cúbica natural")
plot!(p2, xf, Sf, color = :steelblue, lw = 2, label = "S(x)")
scatter!(p2, x0, y0, color = :black, ms = 4, label = "nodos 0:5:100")
p3 = plot(x, max.(err_P, 1e-3), yscale = :log10, color = :firebrick, label = "|P20 - datos|", xlabel = "día",
          ylabel = "error [USD]", title = "(c) Error de representación")
plot!(p3, x, max.(err_S, 1e-3), color = :steelblue, label = "|S - datos|")
fig = plot(p1, p1z, p2, p3, layout = (2, 2), size = (1300, 850), margin = 5Plots.mm)
savefig(fig, joinpath(@__DIR__, "..", "figuras", "P06_splines.png"))
println("Figura guardada en figuras/P06_splines.png")
