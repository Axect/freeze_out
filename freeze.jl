g_eff = 81.
M_pl = 1.221e+19
M = 100.

function rk4(f::Function, x₀::Float64, y₀::Float64, x₁::Float64, n)
    vx = zeros(n + 1)
    vy = zeros(n + 1)
    vx[1] = x = x₀
    vy[1] = y = y₀
    h = (x₁ - x₀) / n
    for i in 1:n
        k₁ = h * f(x, y)
        k₂ = h * f(x + 0.5h, y + 0.5k₁)
        k₃ = h * f(x + 0.5h, y + 0.5k₂)
        k₄ = h * f(x + h, y + k₃)
        vx[i + 1] = x = x₀ + i * h
        vy[i + 1] = y = y + (k₁ + 2k₂ + 2k₃ + k₄) / 6
    end
    return vx, vy
end

function dydx(x::Float64, Y::Float64)
    return -√((π/45.) * g_eff) * M_pl * M / x^2 * 1e-14 * (Y^2 - Y_eq(x)^2)
end

function Y_eq(x::Float64)
    return 45. / (2. * π^4) * √(π / 8.) / g_eff * x^(3. / 2.) * exp(-x)
end

vx, vy = rk4(dydx, 1., Y_eq(1.), 1000., 10000000)

using PyCall

plt = pyimport("pylab")
plt.figure(figsize=(10,6), dpi=300)
plt.loglog(vx, vy)
plt.savefig("plot_julia.png")
