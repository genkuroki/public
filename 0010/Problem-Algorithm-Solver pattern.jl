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
#     display_name: Julia 1.7.0-beta3
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# # Problem-Algorithm-Solver pattern

# %% [markdown]
# In order to use the `(; a, b, c) = p` syntax, require VERSION ≥ v"1.7.0-beta" [#39285](https://github.com/JuliaLang/julia/pull/39285).

# %%
VERSION ≥ v"1.7.0-beta"

# %%
using Plots

# %%
# minimal working example of the problem-algorithm-solver pattern

module FreeFall

"""Problem type"""
Base.@kwdef struct Problem{G, Y0, V0, TS}
    g::G = 9.80665
    y0::Y0 = 0.0
    v0::V0 = 30.0
    tspan::TS = (0.0, 8.0)
end

"""Algorithm types"""
Base.@kwdef struct EulerMethod{T}  dt::T = 0.1 end
Base.@kwdef struct ExactFormula{T} dt::T = 0.1 end
default_algorithm(::Problem) = EulerMethod()

"""Solution type"""
struct Solution{Y, V, T, P<:Problem, A} y::Y; v::V; t::T; prob::P; alg::A end

"""Solver function"""
solve(prob::Problem) = solve(prob, default_algorithm(prob))

function solve(prob::Problem, alg::ExactFormula)
    (; g, y0, v0, tspan) = prob
    (; dt) = alg
    t0, t1 = tspan
    t = range(t0, t1 + dt/2; step = dt)    
    y(t) = y0 + v0*(t - t0) - g*(t - t0)^2/2
    v(t) = v0 - g*(t - t0)
    Solution(y.(t), v.(t), t, prob, alg)
end

function solve(prob::Problem, alg::EulerMethod)
    (; g, y0, v0, tspan) = prob
    (; dt) = alg
    t0, t1 = tspan
    t = range(t0, t1 + dt/2; step = dt)    
    n = length(t)    
    y = Vector{typeof(y0)}(undef, n)
    v = Vector{typeof(v0)}(undef, n)
    y[1] = y0
    v[1] = v0
    for i in 1:n-1
        v[i+1] = v[i] - g*dt
        y[i+1] = y[i] + v[i]*dt
    end
    Solution(y, v, t, prob, alg)
end

end

# %%
earth = FreeFall.Problem()
sol_euler = FreeFall.solve(earth)
sol_exact = FreeFall.solve(earth, FreeFall.ExactFormula())

plot(sol_euler.t, sol_euler.y; label="Euler's method (dt = $(sol_euler.alg.dt))", ls=:auto)
plot!(sol_exact.t, sol_exact.y; label="exact solution", ls=:auto)
title!("On the Earth"; xlabel="t", legend=:bottomleft)

# %%
# change algorithm parameters

earth = FreeFall.Problem()
sol_euler = FreeFall.solve(earth, FreeFall.EulerMethod(dt = 0.01))
sol_exact = FreeFall.solve(earth, FreeFall.ExactFormula())

plot(sol_euler.t, sol_euler.y; label="Euler's method (dt = $(sol_euler.alg.dt))", ls=:auto)
plot!(sol_exact.t, sol_exact.y; label="exact solution", ls=:auto)
title!("On the Earth"; xlabel="t", legend=:bottomleft)

# %%
# change problem parameters

moon = FreeFall.Problem(g = 1.62, tspan = (0.0, 40.0))
sol_euler = FreeFall.solve(moon)
sol_exact = FreeFall.solve(moon, FreeFall.ExactFormula(dt = sol_euler.alg.dt))

plot(sol_euler.t, sol_euler.y; label="Euler's method (dt = $(sol_euler.alg.dt))", ls=:auto)
plot!(sol_exact.t, sol_exact.y; label="exact solution", ls=:auto)
title!("On the Moon"; xlabel="t", legend=:bottomleft)

# %%
# Extensibility of the Julia ecosystem

module SomeExtension

using ..FreeFall
using ..FreeFall: Problem, Solution

Base.@kwdef struct Symplectic2ndOrder{T}  dt::T = 0.1 end

function FreeFall.solve(prob::Problem, alg::Symplectic2ndOrder)
    (; g, y0, v0, tspan) = prob
    (; dt) = alg
    t0, t1 = tspan
    t = range(t0, t1 + dt/2; step = dt)    
    n = length(t)    
    y = Vector{typeof(y0)}(undef, n)
    v = Vector{typeof(v0)}(undef, n)
    y[1] = y0
    v[1] = v0
    for i in 1:n-1
        ytmp = y[i] + v[i]*dt/2
        v[i+1] = v[i] - g*dt
        y[i+1] = ytmp + v[i+1]*dt/2
    end
    Solution(y, v, t, prob, alg)
end

end

earth = FreeFall.Problem()
sol_sympl = FreeFall.solve(earth, SomeExtension.Symplectic2ndOrder(dt = 2.0))
sol_exact = FreeFall.solve(earth, FreeFall.ExactFormula())

plot(sol_sympl.t, sol_sympl.y; label="2nd order symplectic (dt = $(sol_sympl.alg.dt))", ls=:auto)
plot!(sol_exact.t, sol_exact.y; label="exact solution", ls=:auto)
title!("On the Earth"; xlabel="t", legend=:bottomleft)

# %%
# Composability of the Julia ecosystem

using MonteCarloMeasurements

earth = FreeFall.Problem(y0 = 0.0 ± 0.0, v0 = 30.0 ± 1.0)
sol_euler = FreeFall.solve(earth)
sol_sympl = FreeFall.solve(earth, SomeExtension.Symplectic2ndOrder(dt = 2.0))
sol_exact = FreeFall.solve(earth, FreeFall.ExactFormula())

ylim = (-100, 60)
P = plot(sol_euler.t, sol_euler.y; label="Euler's method (dt = $(sol_euler.alg.dt))", ls=:auto)
title!("On the Earth"; xlabel="t", legend=:bottomleft, ylim)
Q = plot(sol_sympl.t, sol_sympl.y; label="2nd order symplectic (dt = $(sol_sympl.alg.dt))", ls=:auto)
title!("On the Earth"; xlabel="t", legend=:bottomleft, ylim)
R = plot(sol_exact.t, sol_exact.y; label="exact solution", ls=:auto)
title!("On the Earth"; xlabel="t", legend=:bottomleft, ylim)
plot(P, Q, R; size=(720, 600))

# %%
