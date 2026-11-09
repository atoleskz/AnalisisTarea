# =============================================================================
#  Problema 4 (Ecuaciones diferenciales 9.31)
#  Problema de valores en la frontera no lineal (barra calentada que irradia)
#  resuelto con el MÉTODO DEL DISPARO + RUNGE-KUTTA DE 4o ORDEN.
#
#     -A k T'' + P σ (T^4 - 273^4) = Q ,   0.5 < x < 1
#     T'(0.5) = 0      (simetría)
#     T(1)    = T_R    (temperatura fija en el extremo)
#
#  Se transforma en el sistema de 1er orden  T' = w,  w' = [Pσ(T^4-273^4) - Q]/(A k)
#  y se ajusta T(0.5) = s hasta que T(1; s) = T_R (método de la secante).
#
#  !!! IMPORTANTE: los valores numéricos de A, k, P y Q deben ser los del
#  !!! ejemplo 9.14 del libro (Nakamura). Aquí están en el bloque CONSTANTES;
#  !!! si los del libro son otros, sólo cámbialos ahí y vuelve a correr.
#
#  Archivo: P04_disparo_rk4.jl         (Tarea de Análisis Numérico 2027-1)
#  Uso:     julia P04_disparo_rk4.jl
#  Salida:  tablas en consola, figuras/P04_disparo.png
# =============================================================================
using Printf
using Plots

# ------------------------------ CONSTANTES -----------------------------------
const A   = 1.0e-4      # área de la sección transversal [m^2]   (barra de 1 cm x 1 cm)
const k   = 60.0        # conductividad térmica [W/(m K)]
const P   = 0.04        # perímetro [m]
const σ   = 5.67e-8     # constante de Stefan-Boltzmann [W/(m^2 K^4)]
const Q   = 50.0        # calor generado por unidad de longitud [W/m]
const Ta  = 273.0       # temperatura del ambiente [K]
const T_R = 273.0       # condición en x = 1 [K]  (extremos a 273 K, como en el ej. 9.14)
const x0, x1 = 0.5, 1.0

# Lado derecho del sistema z = (T, w),  w = T'
G(x, z) = [z[2], (P * σ * (z[1]^4 - Ta^4) - Q) / (A * k)]

# RK4 para sistemas. Devuelve la malla y la solución completa.
function rk4_sistema(G, z0, a, b, h)
    N = round(Int, (b - a) / h)
    xs = [a + i * h for i in 0:N]
    Z = zeros(N + 1, length(z0)); Z[1, :] = z0
    for n in 1:N
        z = Z[n, :]; x = xs[n]
        k1 = h * G(x, z)
        k2 = h * G(x + h / 2, z + k1 / 2)
        k3 = h * G(x + h / 2, z + k2 / 2)
        k4 = h * G(x + h, z + k3)
        Z[n+1, :] = z + (k1 + 2k2 + 2k3 + k4) / 6
    end
    return xs, Z
end

# Residuo del disparo: φ(s) = T(1; s) - T_R
function disparo(s, h)
    xs, Z = rk4_sistema(G, [s, 0.0], x0, x1, h)
    return Z[end, 1] - T_R, xs, Z
end

# --------------------------- 1) Prueba y error -------------------------------
h = 0.01
println("1) Prueba y error (h = $h): se varía T(0.5) = s y se observa T(1)")
@printf("%10s %14s %14s\n", "s = T(0.5)", "T(1; s)", "φ(s)")
for s in (300.0, 350.0, 400.0, 405.0, 410.0, 450.0, 500.0)
    φ, = disparo(s, h)
    @printf("%10.1f %14.4f %14.4f\n", s, φ + T_R, φ)
end

# --------------------------- 2) Bisección + secante --------------------------
# El problema es muy sensible: si T(0.5) es un poco grande, el término T^4 hace
# que la solución del PVI "explote" (T(1) = Inf) antes de llegar a x = 1.
# Por eso primero se acota la raíz con bisección (un disparo que explota cuenta
# como φ > 0) y después se refina con la secante, que converge superlinealmente.
signo(φ) = isfinite(φ) ? sign(φ) : 1.0

function bisec_secante(h; a = 300.0, b = 600.0, tol_b = 1e-3, tol = 1e-9, imprime = true)
    φa, = disparo(a, h)
    imprime && println("\n2) Bisección (h = $h) sobre [$a, $b]")
    imprime && @printf("%4s %16s %16s %16s\n", "it", "a", "b", "φ(m)")
    it = 0
    while b - a > tol_b
        it += 1
        m = (a + b) / 2
        φm, = disparo(m, h)
        imprime && @printf("%4d %16.8f %16.8f %16.6e\n", it, a, b, φm)
        if signo(φm) == signo(φa)
            a, φa = m, φm
        else
            b = m
        end
    end
    # Secante a partir de los extremos del intervalo (ambos con φ finito)
    s0, s1 = a, (a + b) / 2
    φ0, = disparo(s0, h); φ1, = disparo(s1, h)
    if !isfinite(φ1)
        s1 = a + (b - a) / 4; φ1, = disparo(s1, h)
    end
    imprime && println("   Secante:")
    imprime && @printf("%4s %18s %16s\n", "it", "s", "φ(s)")
    for k in 1:30
        s2 = s1 - φ1 * (s1 - s0) / (φ1 - φ0)
        s0, φ0 = s1, φ1
        s1 = s2; φ1, = disparo(s1, h)
        imprime && @printf("%4d %18.10f %16.6e\n", k, s1, φ1)
        abs(φ1) < tol && break
    end
    return s1
end

s_h1 = bisec_secante(0.01)
s_h2 = bisec_secante(0.005, imprime = false)
s_h3 = bisec_secante(0.0025, imprime = false)
@printf("\nT(0.5) con h = 0.01: %.8f K,  h = 0.005: %.8f K,  h = 0.0025: %.8f K\n", s_h1, s_h2, s_h3)
@printf("Orden observado de RK4: log2(|s(h)-s(h/2)| / |s(h/2)-s(h/4)|) = %.2f\n",
        log2(abs(s_h1 - s_h2) / abs(s_h2 - s_h3)))
# Sensibilidad del disparo: dT(1)/ds por diferencia central
δ = 1e-4
sens = (disparo(s_h2 + δ, 0.005)[1] - disparo(s_h2 - δ, 0.005)[1]) / (2δ)
@printf("Sensibilidad dT(1)/dT(0.5) ≈ %.2f  (un error de 0.01 K en T(0.5) cambia T(1) en %.2f K)\n", sens, 0.01sens)

# --------------------------- 3) Perfil y balance -----------------------------
_, xs, Z = disparo(s_h2, 0.005)
println("\n3) Perfil de temperatura (h = 0.005)")
@printf("%8s %14s %14s\n", "x", "T(x) [K]", "T'(x) [K/m]")
for i in 1:10:length(xs)
    @printf("%8.3f %14.6f %14.6f\n", xs[i], Z[i, 1], Z[i, 2])
end

# Balance de energía en la media barra:  Q*0.5 = ∫ Pσ(T^4-273^4) dx  +  (-A k T'(1))
function simpson_datos(y, h)
    n = length(y) - 1
    @assert iseven(n)
    return h / 3 * (y[1] + y[end] + 4sum(y[2:2:end-1]) + 2sum(y[3:2:end-2]))
end
rad = simpson_datos(P * σ .* (Z[:, 1] .^ 4 .- Ta^4), 0.005)
cond = -A * k * Z[end, 2]
@printf("\nBalance: generado Q*0.5 = %.6f W;  radiado = %.6f W;  conducido en x=1 = %.6f W;  suma = %.6f W\n",
        Q * 0.5, rad, cond, rad + cond)

# --------------------------- Gráficas ----------------------------------------
mkpath(joinpath(@__DIR__, "..", "figuras"))
p1 = plot(xlabel = "x [m]", ylabel = "T [K]", title = "Disparos de prueba", legend = :topright)
for s in (round(s_h2) - 10, round(s_h2) - 3, round(s_h2) + 1, round(s_h2) + 3)
    _, xx, ZZ = disparo(s, 0.005)
    Tp = [isfinite(v) && v < 2000 ? v : NaN for v in ZZ[:, 1]]
    plot!(p1, xx, Tp, label = "T(0.5) = $(Int(s)) K", ls = :dash)
end
plot!(p1, xs, Z[:, 1], lw = 3, color = :black, label = @sprintf("solución, T(0.5) = %.2f K", s_h2))
hline!(p1, [T_R], color = :gray, label = "T(1) = $(Int(T_R)) K")
ylims!(p1, (-200, 1000))
# perfil completo usando la simetría respecto a x = 0.5
xfull = vcat(1 .- reverse(xs), xs[2:end]); Tfull = vcat(reverse(Z[:, 1]), Z[2:end, 1])
p2 = plot(xfull, Tfull, lw = 3, color = :firebrick, label = "T(x)", xlabel = "x [m]", ylabel = "T [K]",
          title = "Perfil completo (simetría en x = 0.5)")
vline!(p2, [0.5], ls = :dot, color = :gray, label = "eje de simetría")
ss = range(s_h1 - 30, s_h1 + 4, length = 400)
φs = [disparo(s, 0.01)[1] for s in ss]
φs = [isfinite(v) && abs(v) < 1500 ? v : NaN for v in φs]
p3 = plot(ss, φs, lw = 2, label = "φ(s) = T(1; s) - T_R",
          xlabel = "s = T(0.5) [K]", ylabel = "φ(s) [K]", title = "Función residuo del disparo")
hline!(p3, [0], color = :black, label = "")
scatter!(p3, [s_h1], [0], label = "raíz", color = :firebrick)
fig = plot(p1, p2, p3, layout = (1, 3), size = (1400, 430), margin = 5Plots.mm, left_margin = 10Plots.mm)
savefig(fig, joinpath(@__DIR__, "..", "figuras", "P04_disparo.png"))
println("Figura guardada en figuras/P04_disparo.png")
