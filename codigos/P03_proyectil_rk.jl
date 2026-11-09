# =============================================================================
#  Problema 3 (Ecuaciones diferenciales 9.18)
#  Proyectil con arrastre cuadrático: Runge-Kutta de 2o orden (traducción del
#  programa FORTRAN del libro) y Runge-Kutta de 3er orden.
#
#    u' = -c V u,        u(0) = 150 m/s
#    v' = -g - c V v,    v(0) = 150 m/s,   V = sqrt(u^2 + v^2)
#    x' = u,  y' = v     (x(0) = y(0) = 0),  c = 0.005,  g = 9.8 m/s^2
#
#  Archivo: P03_proyectil_rk.jl        (Tarea de Análisis Numérico 2027-1)
#  Uso:     julia P03_proyectil_rk.jl
#  Salida:  tablas en consola, figuras/P03_*.png
# =============================================================================
using Printf
using Plots

const c = 0.005
const g = 9.8

# --------------------------------------------------------------------------
# a) RK2: traducción línea por línea del programa FORTRAN del libro.
#    Velocidades con RK2 (Heun); posición con regla del trapecio (ecuación B).
# --------------------------------------------------------------------------
function proyectil_rk2(; h = 0.1, nmax = 200, imprime = true)
    UB, VB, T, X, Y = 150.0, 150.0, 0.0, 0.0, 0.0
    hist = [(T, UB, VB, X, Y)]
    imprime && @printf("%12s%12s%12s%12s%12s\n", "TIEMPO", "U", "V", "X", "Y")
    imprime && @printf("%12.6f%12.6f%12.6f%12.6f%12.6f\n", T, UB, VB, X, Y)
    for N in 1:nmax
        T = T + h
        VEL1 = sqrt(UB^2 + VB^2)                       # velocidad absoluta en t(n)
        K1 = -c * VEL1 * UB * h
        L1 = (-g - c * VEL1 * VB) * h
        VEL2 = sqrt((UB + K1)^2 + (VB + L1)^2)         # velocidad absoluta en t(n+1)
        K2 = -c * VEL2 * (UB + K1) * h
        L2 = (-g - c * VEL2 * (VB + L1)) * h
        U = UB + (K1 + K2) / 2
        V = VB + (L1 + L2) / 2
        X = X + 0.5 * (U + UB) * h                     # trapecio: x = ∫ u dt
        Y = Y + 0.5 * (V + VB) * h                     # trapecio: y = ∫ v dt
        UB, VB = U, V
        push!(hist, (T, U, V, X, Y))
        imprime && @printf("%12.6f%12.6f%12.6f%12.6f%12.6f\n", T, U, V, X, Y)
        Y < 0 && break                                 # IF (Y.LT.0) STOP
    end
    return hist
end

# --------------------------------------------------------------------------
# b) RK3 (Kutta, 3er orden) sobre el sistema completo z = (u, v, x, y)
#    k1 = h F(z),  k2 = h F(z + k1/2),  k3 = h F(z - k1 + 2 k2)
#    z_{n+1} = z_n + (k1 + 4 k2 + k3)/6
# --------------------------------------------------------------------------
function F(z)
    u, v = z[1], z[2]
    V = sqrt(u^2 + v^2)
    return [-c * V * u, -g - c * V * v, u, v]
end

function proyectil_rk3(; h = 0.1, nmax = 10^7, imprime = false)
    z = [150.0, 150.0, 0.0, 0.0]
    T = 0.0
    hist = [(T, z...)]
    for N in 1:nmax
        k1 = h * F(z)
        k2 = h * F(z + k1 / 2)
        k3 = h * F(z - k1 + 2k2)
        z = z + (k1 + 4k2 + k3) / 6
        T += h
        push!(hist, (T, z...))
        z[4] < 0 && break
    end
    return hist
end

# Datos de vuelo: altura máxima, tiempo de vuelo y alcance.
# El cruce y = 0 se localiza con interpolación cúbica de Hermite en el último
# paso (usa y' = v y x' = u), para no degradar el orden del método.
function hermite(p0, p1, m0, m1, h, s)          # s ∈ [0,1]
    h00 = 2s^3 - 3s^2 + 1; h10 = s^3 - 2s^2 + s
    h01 = -2s^3 + 3s^2;    h11 = s^3 - s^2
    return h00 * p0 + h10 * h * m0 + h01 * p1 + h11 * h * m1
end
function resumen(hist)
    T = [p[1] for p in hist]; U = [p[2] for p in hist]; V = [p[3] for p in hist]
    X = [p[4] for p in hist]; Y = [p[5] for p in hist]
    imax = argmax(Y)
    n = length(Y); dt = T[n] - T[n-1]
    yH(s) = hermite(Y[n-1], Y[n], V[n-1], V[n], dt, s)
    lo, hi = 0.0, 1.0                              # bisección sobre s
    for _ in 1:60
        mid = (lo + hi) / 2
        yH(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    s = (lo + hi) / 2
    tv = T[n-1] + s * dt
    xr = hermite(X[n-1], X[n], U[n-1], U[n], dt, s)
    return (ymax = Y[imax], tymax = T[imax], tvuelo = tv, alcance = xr)
end

println("a) Programa RK2 (traducción del FORTRAN), h = 0.1  [se imprime cada 1 s]")
h2 = proyectil_rk2(imprime = false)
@printf("%10s%12s%12s%12s%12s\n", "TIEMPO", "U", "V", "X", "Y")
for (i, p) in enumerate(h2)
    (mod(i - 1, 10) == 0 || i == length(h2)) && @printf("%10.2f%12.4f%12.4f%12.4f%12.4f\n", p...)
end
r2 = resumen(h2)

println("\nb) RK3, h = 0.1  [se imprime cada 1 s]")
h3 = proyectil_rk3(h = 0.1)
@printf("%10s%12s%12s%12s%12s\n", "TIEMPO", "U", "V", "X", "Y")
for (i, p) in enumerate(h3)
    (mod(i - 1, 10) == 0 || i == length(h3)) && @printf("%10.2f%12.4f%12.4f%12.4f%12.4f\n", p...)
end
r3 = resumen(h3)

# Referencia: RK3 con h = 1e-4
rref = resumen(proyectil_rk3(h = 1e-4))
# Sin arrastre (c = 0), para comparar: alcance = 2 u0 v0 / g
alc_vacio = 2 * 150 * 150 / g; ymax_vacio = 150^2 / (2g)

println("\nResumen del vuelo")
@printf("%-24s %12s %12s %12s %12s\n", "método", "y_max [m]", "t(y_max)", "t_vuelo [s]", "alcance [m]")
@printf("%-24s %12.4f %12.2f %12.4f %12.4f\n", "RK2 h=0.1", r2.ymax, r2.tymax, r2.tvuelo, r2.alcance)
@printf("%-24s %12.4f %12.2f %12.4f %12.4f\n", "RK3 h=0.1", r3.ymax, r3.tymax, r3.tvuelo, r3.alcance)
@printf("%-24s %12.4f %12.2f %12.4f %12.4f\n", "RK3 h=1e-4 (referencia)", rref.ymax, rref.tymax, rref.tvuelo, rref.alcance)
@printf("%-24s %12.4f %12s %12.4f %12.4f\n", "sin arrastre (c = 0)", ymax_vacio, "-", 2 * 150 / g, alc_vacio)
@printf("\nError en el alcance: RK2 %.3e m,  RK3 %.3e m\n", r2.alcance - rref.alcance, r3.alcance - rref.alcance)

# Convergencia: error máximo en la posición (x, y) sobre los tiempos comunes
# t = 0.4, 0.8, ..., 12.0 (puntos de malla para todos los h considerados)
function estado_en(hist, t)
    i = findmin([abs(p[1] - t) for p in hist])[2]
    return [hist[i][4], hist[i][5]]
end
tcomunes = 0.4:0.4:12.0
ref = proyectil_rk3(h = 1e-4)
errmax(hist) = maximum(maximum(abs.(estado_en(hist, t) .- estado_en(ref, t))) for t in tcomunes)
hs = [0.4, 0.2, 0.1, 0.05, 0.025]
e2 = [errmax(proyectil_rk2(h = h, nmax = 10^6, imprime = false)) for h in hs]
e3 = [errmax(proyectil_rk3(h = h)) for h in hs]
println("\nError máximo en (x, y) para t = 0.4, 0.8, ..., 12 s")
println("   h        error RK2      error RK3   cociente RK2  cociente RK3")
for i in eachindex(hs)
    q2 = i > 1 ? e2[i-1] / e2[i] : NaN; q3 = i > 1 ? e3[i-1] / e3[i] : NaN
    @printf("%6.3f   %12.4e   %12.4e   %10.2f   %10.2f\n", hs[i], e2[i], e3[i], q2, q3)
end
println("Pendientes log-log (orden observado): RK2 ≈ ",
        round((log(e2[1]) - log(e2[end])) / (log(hs[1]) - log(hs[end])), digits = 2),
        ",  RK3 ≈ ", round((log(e3[1]) - log(e3[end])) / (log(hs[1]) - log(hs[end])), digits = 2))

# ------------------------------- Gráficas ------------------------------------
mkpath(joinpath(@__DIR__, "..", "figuras"))
X2 = [p[4] for p in h2]; Y2 = [p[5] for p in h2]
X3 = [p[4] for p in h3]; Y3 = [p[5] for p in h3]
tv = range(0, 2 * 150 / g, length = 300)
p1 = plot(X2, Y2, lw = 2.5, label = "RK2 (programa del libro), h = 0.1", color = :steelblue,
          xlabel = "x [m]", ylabel = "y [m]", title = "9.18 Trayectoria del proyectil", legend = :topright)
plot!(p1, X3, Y3, lw = 2, ls = :dash, label = "RK3, h = 0.1", color = :firebrick)
plot!(p1, 150 .* tv, 150 .* tv .- g / 2 .* tv .^ 2, label = "sin arrastre (c = 0)", color = :gray, ls = :dot)
ylims!(p1, (0, Inf))
T2 = [p[1] for p in h2]
p2 = plot(T2, [p[2] for p in h2], label = "u(t)", lw = 2, xlabel = "t [s]", ylabel = "m/s", title = "Componentes de la velocidad (RK2)")
plot!(p2, T2, [p[3] for p in h2], label = "v(t)", lw = 2)
hline!(p2, [0], color = :black, label = "")
p3 = plot(hs, e2, xscale = :log10, yscale = :log10, marker = :circle, label = "RK2", xlabel = "h", ylabel = "error máx. en (x,y) [m]", title = "Convergencia")
plot!(p3, hs, e3, marker = :square, label = "RK3")
fig = plot(p1, p2, p3, layout = @layout([a{0.55w} [b; c]]), size = (1300, 600), margin = 5Plots.mm)
savefig(fig, joinpath(@__DIR__, "..", "figuras", "P03_proyectil.png"))
println("Figura guardada en figuras/P03_proyectil.png")
