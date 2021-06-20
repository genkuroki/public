# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
@show VERSION
@time using DifferentialEquations
@time using ForwardDiff
@time using Plots

# %%
function f2_mur!(dv, v, u, p, t)
    (; dx) = p
    a, b = firstindex(u), lastindex(u)
    @. @views dv[a+1:b-1] = (u[a:b-2] + u[a+2:b] - 2u[a+1:b-1])/dx^2
    dv[a] = 2(u[a+1] - u[a])/dx^2 - 2v[a]/dx
    dv[b] = 2(u[b-1] - u[b])/dx^2 - 2v[b]/dx
    return
end

# %%
x = range(-10, 10; length=21)
dx = step(x)
p = (; dx)

U(t, x) = 2/3*exp(-(x - t)^2) + 1/3*exp(-(x + t)^2)
V(t, x) = ForwardDiff.derivative(t -> U(t, x), t)
u0 = U.(0, x)
v0 = V.(0, x)
tspan = (0.0, 30.0)

prob = SecondOrderODEProblem(f2_mur!, v0, u0, tspan, p)
sol = solve(prob)

ts = range(sol.prob.tspan...; length=150)
anim = @animate for t in [fill(ts[begin], 10); ts; fill(ts[end], 10)]
    plot(x, sol(t)[end÷2+1:end]; label="", ylim=(-0.25, 1.05), size=(600, 300))
end
gif(anim, "1d_wave_eq_mur_21.gif")

# %%
x = range(-10, 10; length=201)
dx = step(x)
p = (; dx)

U(t, x) = 2/3*exp(-(x - t)^2) + 1/3*exp(-(x + t)^2)
V(t, x) = ForwardDiff.derivative(t -> U(t, x), t)
u0 = U.(0, x)
v0 = V.(0, x)
tspan = (0.0, 20.0)

prob = SecondOrderODEProblem(f2_mur!, v0, u0, tspan, p)
sol = solve(prob)

ts = range(sol.prob.tspan...; length=100)
anim = @animate for t in [fill(ts[begin], 10); ts; fill(ts[end], 10)]
    plot(x, sol(t)[end÷2+1:end]; label="", ylim=(-0.25, 1.05), size=(600, 300))
end
gif(anim, "1d_wave_eq_mur_201.gif")

# %%
