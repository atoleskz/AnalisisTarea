# =============================================================================
#  Problema 2 (Ecuaciones diferenciales 9.1, 9.22, 9.32 y 9.33)
#  Euler hacia adelante (h = 0.5 y h = 0.01), RK4 (h = 0.1),
#  predictor-corrector de Adams de 3er orden (h = 0.5) y de 4o orden (h = 1.0)
#
#  Archivo: P02_euler_rk4_adams.jl     (Tarea de Análisis Numérico 2027-1)
#  Uso:     julia P02_euler_rk4_adams.jl
#  Salida:  tablas en consola, figuras/P02_*.png
# =============================================================================
using Printf
using Plots

# ------------------------------- Problemas -----------------------------------
# f(t, y) para y' = f(t, y)
f_a(t, y) = 1 - t * y                 # a) y' + t y = 1,       y(0) = 1
f_b(t, y) = exp(-t) - 3y              # b) y' + 3y = e^{-t},   y(0) = 1
f_c(t, y) = t^2 - y                   # c) y' = t^2 - y,       y(0) = 0.5
f_d(t, y) = -y * abs(y)               # d) y' + y|y| = 0,      y(0) = 1
f_e(t, y) = sin(t) - sqrt(abs(y))     # e) y' + |y|^{1/2} = sen t, y(0) = 1

# Simpson compuesto (para la solución exacta de a), que contiene una integral)
function simpson(g, a, b, n)
    n = iseven(n) ? n : n + 1
    h = (b - a) / n
    s = g(a) + g(b)
    for i in 1:n-1
        s += (isodd(i) ? 4 : 2) * g(a + i * h)
    end
    return s * h / 3
end

# Soluciones exactas
#  a) y = e^{-t^2/2} (1 + ∫_0^t e^{s^2/2} ds)   (la integral se evalúa con Simpson, n = 20000)
#  b) y = (e^{-t} + e^{-3t})/2
#  c) y = t^2 - 2t + 2 - 1.5 e^{-t}
#  d) y = 1/(1 + t)
#  e) sin forma cerrada: referencia con RK4 y h = 1e-5
y_a(t) = exp(-t^2 / 2) * (1 + (t == 0 ? 0.0 : simpson(s -> exp(s^2 / 2), 0.0, t, 20000)))
y_b(t) = (exp(-t) + exp(-3t)) / 2
y_c(t) = t^2 - 2t + 2 - 1.5exp(-t)
y_d(t) = 1 / (1 + t)

# ------------------------------- Métodos -------------------------------------
function euler(f, y0, t0, tf, h)
    N = round(Int, (tf - t0) / h)
    t = [t0 + i * h for i in 0:N]
    y = zeros(N + 1); y[1] = y0
    for n in 1:N
        y[n+1] = y[n] + h * f(t[n], y[n])
    end
    return t, y
end

function rk4(f, y0, t0, tf, h)
    N = round(Int, (tf - t0) / h)
    t = [t0 + i * h for i in 0:N]
    y = zeros(N + 1); y[1] = y0
    for n in 1:N
        k1 = h * f(t[n], y[n])
        k2 = h * f(t[n] + h / 2, y[n] + k1 / 2)
        k3 = h * f(t[n] + h / 2, y[n] + k2 / 2)
        k4 = h * f(t[n] + h, y[n] + k3)
        y[n+1] = y[n] + (k1 + 2k2 + 2k3 + k4) / 6
    end
    return t, y
end

# Predictor-corrector de Adams. Arranque con RK4 (como sugiere el problema 9.36).
#  Orden 3:  P (Adams-Bashforth): y_{n+1} = y_n + h/12 (23 f_n - 16 f_{n-1} + 5 f_{n-2})
#            C (Adams-Moulton):   y_{n+1} = y_n + h/12 (5 f_{n+1} + 8 f_n - f_{n-1})
#  Orden 4:  P: y_{n+1} = y_n + h/24 (55 f_n - 59 f_{n-1} + 37 f_{n-2} - 9 f_{n-3})
#            C: y_{n+1} = y_n + h/24 (9 f_{n+1} + 19 f_n - 5 f_{n-1} + f_{n-2})
function adams_pc(f, y0, t0, tf, h, orden)
    N = round(Int, (tf - t0) / h)
    t = [t0 + i * h for i in 0:N]
    y = zeros(N + 1); y[1] = y0
    m = orden - 1                                  # pasos de arranque
    _, yr = rk4(f, y0, t0, t0 + min(m, N) * h, h)
    y[1:min(m, N)+1] .= yr
    F = [f(t[i], y[i]) for i in 1:min(m, N)+1]
    for n in (m+1):N                               # calcula y[n+1]
        if orden == 3
            yp = y[n] + h / 12 * (23F[n] - 16F[n-1] + 5F[n-2])
            fp = f(t[n+1], yp)
            y[n+1] = y[n] + h / 12 * (5fp + 8F[n] - F[n-1])
        else
            yp = y[n] + h / 24 * (55F[n] - 59F[n-1] + 37F[n-2] - 9F[n-3])
            fp = f(t[n+1], yp)
            y[n+1] = y[n] + h / 24 * (9fp + 19F[n] - 5F[n-1] + F[n-2])
        end
        push!(F, f(t[n+1], y[n+1]))
    end
    return t, y
end

# ------------------------------- Ejecución -----------------------------------
# Referencia para e): RK4 con h = 1e-5, interpolada en la malla pedida
tref_e, yref_e = rk4(f_e, 1.0, 0.0, 5.0, 1e-5)
y_e(t) = yref_e[round(Int, t / 1e-5) + 1]

problemas = [("a", f_a, 1.0, y_a), ("b", f_b, 1.0, y_b), ("c", f_c, 0.5, y_c),
             ("d", f_d, 1.0, y_d), ("e", f_e, 1.0, y_e)]

valor_en(t, y, ts) = y[findmin(abs.(t .- ts))[2]]
tabla_t = 0.0:1.0:5.0

println("="^78)
println(" 9.1  Euler hacia adelante, h = 0.5 (tabla completa, la que se hace a mano)")
println("="^78)
for (nom, f, y0, yex) in problemas
    t, y = euler(f, y0, 0.0, 5.0, 0.5)
    println("\nInciso $nom)")
    @printf("%6s %14s %14s %12s\n", "t", "Euler h=0.5", "exacta", "error")
    for i in eachindex(t)
        @printf("%6.2f %14.6f %14.6f %12.4e\n", t[i], y[i], yex(t[i]), y[i] - yex(t[i]))
    end
end

println("\n", "="^78)
println(" Resumen en t = 1,...,5: errores de cada método (valor numérico - exacto)")
println("="^78)
resultados = Dict{String,Any}()
for (nom, f, y0, yex) in problemas
    sol = Dict(
        "Euler h=0.5"  => euler(f, y0, 0.0, 5.0, 0.5),
        "Euler h=0.01" => euler(f, y0, 0.0, 5.0, 0.01),
        "RK4 h=0.1"    => rk4(f, y0, 0.0, 5.0, 0.1),
    )
    if nom in ("a", "b", "c")
        sol["Adams3 h=0.5"] = adams_pc(f, y0, 0.0, 5.0, 0.5, 3)
        sol["Adams4 h=1.0"] = adams_pc(f, y0, 0.0, 5.0, 1.0, 4)
    end
    resultados[nom] = sol
    metodos = filter(k -> haskey(sol, k), ["Euler h=0.5", "Euler h=0.01", "RK4 h=0.1", "Adams3 h=0.5", "Adams4 h=1.0"])
    println("\nInciso $nom)")
    @printf("%4s %11s", "t", "exacta")
    for m in metodos; @printf(" %13s", m); end
    println()
    for ts in tabla_t
        @printf("%4.1f %11.6f", ts, yex(ts))
        for m in metodos
            t, y = sol[m]
            @printf(" %13.3e", valor_en(t, y, ts) - yex(ts))
        end
        println()
    end
    # Error máximo en la malla de cada método
    for m in metodos
        t, y = sol[m]
        emax = maximum(abs.(y .- yex.(t)))
        @printf("   error máximo %-13s = %.4e\n", m, emax)
    end
end

# ------------------------------- Gráficas ------------------------------------
mkpath(joinpath(@__DIR__, "..", "figuras"))
tt = range(0, 5, length = 501)
plots_sol = []
for (nom, f, y0, yex) in problemas
    sol = resultados[nom]
    p = plot(tt, yex.(tt), lw = 2.5, color = :black, label = "exacta", title = "9.1 $nom)", xlabel = "t", ylabel = "y")
    t, y = sol["Euler h=0.5"];  plot!(p, t, y, marker = :circle, ms = 3, label = "Euler h=0.5", color = :firebrick)
    t, y = sol["Euler h=0.01"]; plot!(p, t, y, label = "Euler h=0.01", color = :steelblue, ls = :dash)
    t, y = sol["RK4 h=0.1"];    plot!(p, t, y, label = "RK4 h=0.1", color = :seagreen, ls = :dot, lw = 2)
    if haskey(sol, "Adams3 h=0.5")
        t, y = sol["Adams3 h=0.5"]; plot!(p, t, y, marker = :diamond, ms = 3, label = "Adams3 h=0.5", color = :purple)
        t, y = sol["Adams4 h=1.0"]; plot!(p, t, y, marker = :square, ms = 3, label = "Adams4 h=1", color = :darkorange)
    end
    push!(plots_sol, p)
end
fig1 = plot(plots_sol..., layout = (2, 3), size = (1300, 750), legend = :best, legendfontsize = 6, margin = 4Plots.mm)
savefig(fig1, joinpath(@__DIR__, "..", "figuras", "P02_soluciones.png"))

plots_err = []
for (nom, f, y0, yex) in problemas
    sol = resultados[nom]
    p = plot(title = "Error |y_n - y(t_n)|, inciso $nom)", xlabel = "t", yscale = :log10, legendfontsize = 6)
    for (m, c) in [("Euler h=0.5", :firebrick), ("Euler h=0.01", :steelblue), ("RK4 h=0.1", :seagreen),
                   ("Adams3 h=0.5", :purple), ("Adams4 h=1.0", :darkorange)]
        haskey(sol, m) || continue
        t, y = sol[m]
        e = max.(abs.(y .- yex.(t)), 1e-16)
        plot!(p, t[2:end], e[2:end], label = m, color = c, marker = length(t) < 20 ? :circle : :none, ms = 3)
    end
    push!(plots_err, p)
end
fig2 = plot(plots_err..., layout = (2, 3), size = (1300, 750), margin = 4Plots.mm)
savefig(fig2, joinpath(@__DIR__, "..", "figuras", "P02_errores.png"))
println("\nFiguras guardadas en figuras/P02_soluciones.png y figuras/P02_errores.png")
