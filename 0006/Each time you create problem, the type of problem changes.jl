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

# %% [markdown]
# https://discourse.julialang.org/t/function-compiles-every-time/63273

# %%
using StaticArrays
using AstrodynamicalModels
using DifferentialEquations

# %%
using StaticArrays
using AstrodynamicalModels
using DifferentialEquations

defaults = (; reltol=1e-14, abstol=1e-14, save_everystep=false)
kwargs   = (;)
options  = merge(defaults, kwargs)

r = randn(3) * 1e4
v = randn(3) * 1e3
T = 10_000.0

p1 = ODEProblem(R2BP, MVector{6}(vcat(r,v)), (0.0, T), MVector{1}(398600.435))
T1 = typeof(p1)

p2 = ODEProblem(R2BP, MVector{6}(vcat(r,v)), (0.0, T), MVector{1}(398600.435))
T2 = typeof(p2)

T1 === T2 #> false

# %%
using StaticArrays
using AstrodynamicalModels
using DifferentialEquations

function propagate(r, v, T; kwargs...)

  defaults = (; reltol=1e-14, abstol=1e-14, save_everystep=false)
  options  = merge(defaults, kwargs)

  problem   = ODEProblem(R2BP, MVector{6}(vcat(r,v)), (0.0, T), MVector{1}(398600.435))
  solutions = solve(problem; options...)

end

r = randn(3) * 1e4
v = randn(3) * 1e3
T = 10_000.0

@time propagate(r,v,T) #  28.067676 seconds (93.98 M allocations: 5.404 GiB, 3.61% gc time)
@time propagate(r,v,T) #   7.187065 seconds (13.05 M allocations: 704.718 MiB, 1.34% gc time, 99.70% compilation time)
@time propagate(r,v,T) #   7.154281 seconds (13.06 M allocations: 705.808 MiB, 1.79% gc time, 99.43% compilation time)

# %%
using StaticArrays
using AstrodynamicalModels
using DifferentialEquations

defaults = (; reltol=1e-14, abstol=1e-14, save_everystep=false)
kwargs   = (;)
options  = merge(defaults, kwargs)

r = randn(3) * 1e4
v = randn(3) * 1e3
T = 10_000.0

problem   = ODEProblem(R2BP, MVector{6}(vcat(r,v)), (0.0, T), MVector{1}(398600.435))

@time solutions = solve(problem; options...) #  15.728376 seconds (50.96 M allocations: 2.763 GiB, 4.17% gc time, 2.34% compilation time)
@time solutions = solve(problem; options...) #   0.012221 seconds (540.66 k allocations: 9.171 MiB)
@time solutions = solve(problem; options...) #   0.012453 seconds (540.66 k allocations: 9.171 MiB)

# %%
using StaticArrays
using AstrodynamicalModels
using DifferentialEquations

function propagate_remake(problem_default, r, v, T; kwargs...)

  defaults = (; reltol=1e-14, abstol=1e-14, save_everystep=false)
  options  = merge(defaults, kwargs)

  problem_remaked = remake(problem; 
        u0 = MVector{6}(vcat(r,v)), tspan = (0.0, T), p = MVector{1}(398600.435))
  solutions = solve(problem_remaked; options...)

end

r = randn(3) * 1e4
v = randn(3) * 1e3
T = 10_000.0
problem = ODEProblem(R2BP, MVector{6}(vcat(r,v)), (0.0, T), MVector{1}(398600.435))

r = randn(3) * 1e4
v = randn(3) * 1e3
T = 10_000.0
@time sol1 = propagate_remake(problem, r, v, T)

r = randn(3) * 1e4
v = randn(3) * 1e3
T = 10_000.0
@time sol2 = propagate_remake(problem, r, v, T)

r = randn(3) * 1e4
v = randn(3) * 1e3
T = 10_000.0
@time sol3 = propagate_remake(problem, r, v, T);

# 9.130384 seconds (13.13 M allocations: 715.117 MiB, 1.53% gc time, 0.16% compilation time)
# 0.014155 seconds (452.63 k allocations: 7.765 MiB)
# 0.013558 seconds (448.20 k allocations: 7.689 MiB)

# %%
using Plots
PP = []
for sol in (sol1, sol2, sol3)
    t = range(sol.prob.tspan...; length=1000)
    y = hcat((sol.(t; idxs=i) for i in 1:3)...)
    P = plot(t, y; label="")
    push!(PP, P)
end
plot(PP...; size=(720, 150), layout=(1, 3))

# %%
