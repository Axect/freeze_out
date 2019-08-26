g_eff = 81.
M_pl = 1.221e+19
M = 100.

function rk4(f::Function, x₀::Float64, y₀::Float64, x₁::Float64, sv, n)
    vx = zeros(n + 1)
    vy = zeros(n + 1)
    vx[1] = x = x₀
    vy[1] = y = y₀
    h = (x₁ - x₀) / n
    for i in 1:n
        k₁ = h * f(x, y, sv)
        k₂ = h * f(x + 0.5h, y + 0.5k₁, sv)
        k₃ = h * f(x + 0.5h, y + 0.5k₂, sv)
        k₄ = h * f(x + h, y + k₃, sv)
        vx[i + 1] = x = x₀ + i * h
        vy[i + 1] = y = y + (k₁ + 2k₂ + 2k₃ + k₄) / 6
    end
    return vx, vy
end

function dydx(x::Float64, Y::Float64, sv)
    return -√((π/45.) * g_eff) * M_pl * M / x^2 * sv * (Y^2 - Y_eq(x)^2)
end

function Y_eq(x::Float64)
    return 45. / (2. * π^4) * √(π / 8.) / g_eff * x^(3. / 2.) * exp(-x)
end

vx, vy = rk4(dydx, 1., Y_eq(1.), 100., 1e-14, 10000000)
#vx2, vy2 = rk4(dydx, 1., Y_eq(1.), 100., 1e-15, 10000000)
vx3, vy3 = rk4(dydx, 1., Y_eq(1.), 100., 1e-16, 10000000)
vy_eq = map(Y_eq, vx)

using PyCall

plt = pyimport("pylab")

plt.rc("text", usetex=true)
plt.rc("font", family="serif")

plt.figure(figsize=(10,6), dpi=300)
plt.title("relic abundance")
plt.xlabel("\$x=M/T\$")
plt.ylabel("\$Y\$")
plt.loglog(vx, vy)
#plt.loglog(vx, vy2)
plt.loglog(vx, vy3)
plt.loglog(vx, vy_eq)

axes = plt.gca()
axes.set_ylim([1e-7, 1e-3])

plt.savefig("plot_julia.png")
