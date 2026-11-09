# =============================================================================
#  Problema 1 (Problemas de errores, ejercicio 11)
#  Esquema de Horner: análisis del error de redondeo y cota a posteriori
#
#  Archivo: P01_horner_error.jl        (Tarea de Análisis Numérico 2027-1)
#  Uso:     julia P01_horner_error.jl
#  Salida:  tabla en consola y figuras/P01_horner_error.png
# =============================================================================
using Printf
using Plots

# -----------------------------------------------------------------------------
# Horner con cota de error dinámica (recursión e_r del ejercicio 11(b)).
#   b_n := a_n,  b_r := x*b_{r+1} + a_r   (r = n-1,...,0)
#   e_n := |a_n|/2,  e_r := |x|*e_{r+1} + |b'_r|
#   |A(x) - b'_0| <= eps/(1-eps) * (2 e_0 - |b'_0|)
# a = [a_0, a_1, ..., a_n]  (Julia indexa desde 1: a[k+1] = a_k)
# -----------------------------------------------------------------------------
function horner_con_cota(a::Vector{T}, x::T) where {T<:AbstractFloat}
    n = length(a) - 1
    u = eps(T) / 2                 # unidad de redondeo (eps del problema)
    b = a[n+1]                     # b'_n
    e = abs(a[n+1]) / 2            # e_n
    for r in n-1:-1:0
        b = x * b + a[r+1]         # b'_r = fl(fl(x b'_{r+1}) + a_r)
        e = abs(x) * e + abs(b)    # e_r
    end
    cota = u / (1 - u) * (2e - abs(b))
    return b, cota
end

# Horner "exacto" con aritmética de precisión extendida (referencia)
function horner_exacto(a::Vector, x)
    ab = BigFloat.(a); xb = BigFloat(x)
    b = ab[end]
    for r in length(ab)-1:-1:1
        b = xb * b + ab[r]
    end
    return b
end

# Verificación numérica del inciso (a): A(z) = (z - x) B(z) + b_0
function verifica_a(a, x, z)
    n = length(a) - 1
    b = zeros(BigFloat, n + 1)          # b[r+1] = b_r
    b[n+1] = a[n+1]
    for r in n-1:-1:0
        b[r+1] = x * b[r+2] + a[r+1]
    end
    A = sum(a[r+1] * big(z)^r for r in 0:n)
    B = sum(b[r+1] * big(z)^(r-1) for r in 1:n)
    return A, (z - x) * B + b[1]
end

# -----------------------------------------------------------------------------
# Polinomio de prueba: A(x) = (x - 2)^9 expandido (cancelación fuerte cerca de 2)
# -----------------------------------------------------------------------------
binom(n, k) = factorial(big(n)) ÷ (factorial(big(k)) * factorial(big(n - k)))
grado = 9
a = Float64[binom(grado, k) * (-2.0)^(grado - k) for k in 0:grado]   # a_k

println("Inciso (a): verificación de A(z) = (z - x)B(z) + b0 con x = 1.5")
for z in (-1.0, 0.3, 2.7)
    A, D = verifica_a(BigFloat.(a), big(1.5), z)
    @printf("  z = %5.2f   A(z) = %+.10e   (z-x)B(z)+b0 = %+.10e\n", z, Float64(A), Float64(D))
end

println("\nInciso (b): error real vs. cota  |A(x) - b'_0| <= eps/(1-eps)(2e_0 - |b'_0|)")
@printf("%8s %15s %15s %13s %13s %6s\n", "x", "b'_0 (Float64)", "A(x) exacto", "error real", "cota", "¿ok?")
xs_tabla = [1.90, 1.95, 1.99, 2.00, 2.01, 2.05, 2.10]
for x in xs_tabla
    b0, cota = horner_con_cota(a, x)
    ex = horner_exacto(a, x)
    err = abs(Float64(ex - big(b0)))
    @printf("%8.3f %15.6e %15.6e %13.4e %13.4e %6s\n", x, b0, Float64(ex), err, cota, err <= cota ? "sí" : "NO")
end

# Barrido fino para la gráfica
xs = range(1.92, 2.08, length = 801)
errs = Float64[]; cotas = Float64[]; vals = Float64[]
for x in xs
    b0, c = horner_con_cota(a, x)
    push!(vals, b0)
    push!(errs, max(abs(Float64(horner_exacto(a, x) - big(b0))), 1e-30))
    push!(cotas, c)
end
@printf("\nMáximo error real en [1.92, 2.08]: %.3e   máxima cota: %.3e\n", maximum(errs), maximum(cotas))
@printf("Proporción de puntos donde error <= cota: %d / %d\n", count(errs .<= cotas), length(xs))

# Cota a priori (clásica):  |A(x)-b'_0| <= 2n*u * sum |a_k||x|^k  (aprox.)
u = eps(Float64) / 2
apriori = [2grado * u * sum(abs(a[k+1]) * abs(x)^k for k in 0:grado) for x in xs]

p1 = plot(xs, vals, label = "Horner en Float64", lw = 1, color = :steelblue,
          xlabel = "x", ylabel = "A(x)", title = "A(x) = (x-2)^9 expandido")
plot!(p1, xs, (xs .- 2) .^ 9, label = "(x-2)^9 exacto", lw = 2, color = :black, ls = :dash)
p2 = plot(xs, errs, yscale = :log10, label = "error real", color = :steelblue,
          xlabel = "x", ylabel = "error absoluto", title = "Error vs. cotas", legend = :bottomright)
plot!(p2, xs, cotas, label = "cota a posteriori (ej. 11)", lw = 2, color = :firebrick)
plot!(p2, xs, apriori, label = "cota a priori", lw = 2, color = :darkorange, ls = :dash)
fig = plot(p1, p2, layout = (1, 2), size = (1100, 420), margin = 5Plots.mm)
mkpath(joinpath(@__DIR__, "..", "figuras"))
savefig(fig, joinpath(@__DIR__, "..", "figuras", "P01_horner_error.png"))
println("Figura guardada en figuras/P01_horner_error.png")
