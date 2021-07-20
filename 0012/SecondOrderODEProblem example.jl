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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://discourse.julialang.org/t/empty-container-reduce-error-diff-eq/64956

# %%
using DifferentialEquations
using LinearAlgebra
using Plots

# %%
function plot_2ndorder(sol; l1=:bottomright, l2=:bottomright, size=(720, 300), kwargs...)
    t = range(sol.prob.tspan...; length=400)
    n = Base.size(sol, 1) ÷ 2
    v = vcat((t -> sol(t)[1:n]').(t)...)
    d = vcat((t -> sol(t)[n+1:end]').(t)...)
    P = plot(t, d; label=permutedims(["\$d_{$i}\$" for i in 1:n]), legend=l1)
    Q = plot(t, v; label=permutedims(["\$v_{$i}\$" for i in 1:n]), legend=l2)
    plot(P, Q; size, kwargs...)
end

# %%
function g_naive!(dv, v, d, p, t)
    m, f, c, K = p
    # The following is equivalent to dv .= M \ (f - (C * v + K * d)),
    # where K = Matrix, C = Diagonal(c), M = Diagonal(m)
    dv .= m .\ (f .- (c .* v .+ K * d))
    return
end

# %%
m = Float64[1, 2, 3, 4]
f = Float64[-10, -10, -10, -10]
c = Float64[1, 1, 1, 1]
K = Float64[
     2 -1  0  0
    -1  2 -1  0
     0  1  2 -1
     0  0 -1  2
]
p = (m, f, c, K)
v0 = Float64[0, 0, 0, 0]
d0 = Float64[3, 1, -1, -3]
tspan = (0.0, 10.0)

# %%
prob = SecondOrderODEProblem(g_naive!, v0, d0, tspan, p)
sol = solve(prob)
plot_2ndorder(sol)

# %%
function g!(dv, v, d, p, t)
    m, f, c, K = p
    # The following is equivalent to dv .= M \ (f - (C * v + K * d)),
    # where K = Matrix, C = Diagonal(c), M = Diagonal(m)
    mul!(dv, K, d)
    @. dv = m \ (f - (c * v + dv))
    return
end

# %%
prob = SecondOrderODEProblem(g!, v0, d0, tspan, p)
sol = solve(prob)
plot_2ndorder(sol)

# %%
using BenchmarkTools

@btime g_naive!($(similar(v0)), $v0, $d0, $p, 0.0)
@btime g!($(similar(v0)), $v0, $d0, $p, 0.0)

# %%
