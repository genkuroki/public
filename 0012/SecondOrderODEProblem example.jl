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
K = [
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
@time sol = solve(prob)
@time sol = solve(prob)
@time sol = solve(prob)

t = range(sol.prob.tspan...; length=400)
v = vcat((t -> sol(t)[1:4]').(t)...)
d = vcat((t -> sol(t)[5:8]').(t)...)
P = plot(t, d; label=permutedims(["\$d_$i\$" for i in 1:4]), legend=:bottomleft)
Q = plot(t, v; label=permutedims(["\$v_$i\$" for i in 1:4]), legend=:bottomright)
plot(P, Q; size=(720, 300))

# %%
function g!(dv, v, d, p, t)
    m, f, c, K = p
    # The following is equivalent to dv .= M \ (f - (C * v + K * d)),
    # where K = Matrix, C = Diagonal(c), M = Diagonal(m)
    mul!(dv, K, d)
    @. dv += c * v
    @. dv = f - dv
    @. dv = m \ dv
    return
end

# %%
prob = SecondOrderODEProblem(g!, v0, d0, tspan, p)
@time sol = solve(prob)
@time sol = solve(prob)
@time sol = solve(prob)

t = range(sol.prob.tspan...; length=400)
v = vcat((t -> sol(t)[1:4]').(t)...)
d = vcat((t -> sol(t)[5:8]').(t)...)
P = plot(t, d; label=permutedims(["\$d_$i\$" for i in 1:4]), legend=:bottomleft)
Q = plot(t, v; label=permutedims(["\$v_$i\$" for i in 1:4]), legend=:bottomright)
plot(P, Q; size=(720, 300))

# %%
