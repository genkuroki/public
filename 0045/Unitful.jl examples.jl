# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Unitful
cm = u"cm"
mm = u"mm"

4cm ÷ 5mm

# %%
4cm / 5mm

# %%
4cm / 5mm |> NoUnits

# %%
4cm / 5mm + 0

# %%
using Unitful
m, kg, s = u"m", u"kg", u"s"
using Plots
default(fmt=:png)

"""du/dt = F(u, v), dv/dt = G(u, v)"""
function symplectic_euler(G, F, v0, u0, tspan, param, dt)
    t = range(tspan...; step = dt)
    u = similar(t, typeof(u0))
    v = similar(t, typeof(v0))
    u[begin], v[begin] = u0, v0
    for i in eachindex(t)[begin:end-1]
        v[i+1] = v[i] + G(u[i], v[i], param)*dt
        u[i+1] = u[i] + F(u[i], v[i+1], param)*dt
    end
    t, u, v
end

F(θ, p, (; M, L)) = p/(M*L^2)
G(θ, p, (; M, L, g)) = -M*g*L*sin(θ)
param = (M = 0.5kg, L = 2.0m, g = 9.8m/s^2)
ymax = 0.99
θ0 = 0.0
p0 = ((; M, L, g) = param; M*L^2*√(g/L)*2ymax)
tspan = (0.0s, 20.0s)
dt = 0.01s

t, θ, p = symplectic_euler(G, F, p0, θ0, tspan, param, dt)

P1 = plot(t, θ; label="", xguide="t", yguide="θ", c=1)
P2 = plot(t, p; label="", xguide="t", yguide="p", c=2)
plot(P1, P2; size=(600, 500), layout=(2, 1))

# %%
Any[t θ p]

# %%
using Elliptic: Jacobi
sol_exact(t, (; L, g), ymax) = 2asin(ymax*Jacobi.sn(√(g/L)*t, ymax^2))

plot(t, θ; label="numerical", xguide="t", yguide="θ")
plot!(t, t -> sol_exact(t, param, ymax); label="exact", ls=:dash)
plot!(size=(600, 250), legend=:topright)
plot!(leftmargin=4Plots.mm, bottommargin=4Plots.mm)

# %%
