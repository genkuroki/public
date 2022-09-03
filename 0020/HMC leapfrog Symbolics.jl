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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# # Hamiltonian Monte Carlo with leapfrog
#
# * 黒木玄
# * 2021-12-06, 2022-09-03
#
# Scalar version: https://github.com/genkuroki/public/blob/main/0018/HMC%20leapfrog.ipynb
#
# 2022-09-03: カーネルをJulia v1.8.0に変更して, ConcreteStructs.jl と Parameters.jl への依存を無くした.

# %% tags=[]
# Sample code of Hamiltonian Monte Carlo with leapfrog

module My

using LinearAlgebra: dot
using ForwardDiff: gradient
using Random: default_rng, randn!
using StaticArrays: SVector, MVector

# Type of Leapfrog Problem
struct LFProblem{dim, F1, F2, F3, T, I}
    ϕ::F1
    H::F2
    F::F3
    dt::T
    nsteps::I
end

"""
    LFProblem(dim, ϕ, H, F, dt, nsteps)

Assume that `ϕ(x, param)` is a potential function of a `dim`-dimensional vector `x` and a parameter `param`, `H(x, v, param)` and `F(x, param)` the Hamiltonian and force functions corresponding to `ϕ`, `dt` a step of discretized time, and `nsteps` the number of steps.  Then it returns the Problem object of solving the Hamiltonian equation with leapfrog method.
"""
function LFProblem(dim, ϕ, H, F, dt, nsteps)
    F1, F2, F3, T, I = typeof(ϕ), typeof(H), typeof(F), typeof(dt), typeof(nsteps)
    LFProblem{dim, F1, F2, F3, T, I}(ϕ, H, F, dt, nsteps)
end

"""
    LFProblem(dim, ϕ; dt = 1.0, nsteps = 40)

Assume that `ϕ(x, param)` is a potential function of a `dim`-dimensional vector `x` and a parameter `param`.  Then it defines the Hamiltonian function `H(x, v, param)` and the force function `F(x, param)` corresponding to `ϕ` and returns `LFProblem(dim, ϕ, H, F, dt, nsteps)`.
"""
function LFProblem(dim, ϕ; dt = 1.0, nsteps = 40)
    H(x, v, param) = dot(v, v)/2 + ϕ(x, param)
    F(x, param) = -gradient(x -> ϕ(x, param), x)
    LFProblem(dim, ϕ, H, F, dt, nsteps)
end

"""
    LFProblem(dim, ϕ, ∇ϕ; dt = 1.0, nsteps = 40)

Assume that `ϕ(x, param)` is a potential function of a `dim`-dimensional vector `x` and a parameter `param` and `∇ϕ` its gradient with respect to `x`.  Then it deifnes the Hamiltonian function `H(x, v, param)` and the force function `F(x, param)` corresponding to `ϕ` and `∇ϕ` and returns `LFProblem(dim, ϕ, H, F, dt, nsteps)`.
"""
function LFProblem(dim, ϕ, ∇ϕ; dt = 1.0, nsteps = 40)
    H(x, v, param) = dot(v, v)/2 + ϕ(x, param)
    F(x, param) = -∇ϕ(x, param)
    LFProblem(dim, ϕ, H, F, dt, nsteps)
end

"""
    solve(lf::LFProblem, x, v, param)

numerically solve the Hamilton's equation of motion given by `lf` with leapfrog method, where (`x`, `v`) is the initial value and `param` is the parameter of the potential function `lf.ϕ`.
"""
function solve(lf::LFProblem, x, v, param)
    (; F, dt, nsteps) = lf
    v = v + F(x, param)*dt/2
    x = x + v*dt
    for _ in 2:nsteps
        v = v + F(x, param)*dt
        x = x + v*dt
    end
    v = v + F(x, param)*dt/2
    x, v
end

@inline function _update!(lf::LFProblem{dim}, x, vtmp, param, rng) where dim
    (; H) = lf
    v = SVector{dim}(randn!(rng, vtmp))
    xnew, vnew = solve(lf, x, v, param)
    dH = H(xnew, vnew, param) - H(x, v, param)
    rand(rng) ≤ exp(-dH) ? xnew : x
end

"""
    HMC(lf::LFProblem{dim}, param = nothing;
        niters = 10^5, thin = 1, nwarmups = 0, rng = default_rng(),
        init = SVector{dim}(randn(rng, dim))) where dim

generates the sample of the distribution given by the probability density function proportioal to exp(-`lf.ϕ(x, param)`) by Hamiltonian Monte Carlo method.
"""
function HMC(lf::LFProblem{dim}, param = nothing;
        niters = 10^5, thin = 1, nwarmups = 0, rng = default_rng(),
        init = SVector{dim}(randn(rng, dim))) where dim
    vtmp = MVector{dim}(zeros(eltype(init), dim))
    x = init
    for _ in 1:nwarmups
        x = _update!(lf, x, vtmp, param, rng)
    end
    sample = Vector{typeof(init)}(undef, niters)
    for i in 1:niters
        for _ in 1:thin
            x = _update!(lf, x, vtmp, param, rng)
        end
        @inbounds sample[i] = x
    end
    sample
end

end

# %%
?My.LFProblem

# %%
?My.solve

# %%
?My.HMC

# %%
using Plots
using BenchmarkTools
using StaticArrays
using LinearAlgebra
using KernelDensity
using Statistics
using QuadGK
using Distributions
using Symbolics

# %% [markdown]
# ## 2-dimensional normal distribution

# %%
A = @SMatrix [
     1  1/2
    1/2  1
]
param = (; A = A)
ϕ(x, param) = dot(x, param.A, x)/2
lf = My.LFProblem(2, ϕ)

# %%
@time sample = My.HMC(lf, param)
@time sample = My.HMC(lf, param)
@time sample = My.HMC(lf, param);

# %%
@btime My.HMC($lf, $param);

# %%
X, Y = first.(sample), last.(sample)
d = InterpKDE(kde((X, Y)))
x, y = range(extrema(X)...; length=201), range(extrema(Y)...; length=201)
heatmap(x, y, (x, y) -> pdf(d, x, y); size=(450, 400), right_margin=3Plots.mm)

# %%
f(n) = mean(x -> x*x', @view sample[1:n])

n = 1:1000
S = f.(n)
S11 = (S -> S[1,1]).(S)
S22 = (S -> S[2,2]).(S)
S12 = (S -> S[1,2]).(S)

ymin = min(-1.5, minimum(S11), minimum(S22), minimum(S12))
ymax = max(2.5, maximum(S11), maximum(S22), maximum(S12))

plot(ylim = (ymin, ymax))
plot!(S11; label="s11", c=1)
hline!([inv(A)[1,1]]; label="", c=1, ls=:dash)
plot!(S22; label="s22", c=2)
hline!([inv(A)[2,2]]; label="", c=2, ls=:dash)
plot!(S12; label="s12", c=3)
hline!([inv(A)[1,2]]; label="", c=3, ls=:dash)

# %% [markdown]
# ## φ(x) = a(x² - 1)²

# %%
ϕ4(x, a) = a * (x[1]^2 - 1)^2
a = [3, 4, 5, 6, 7, 8]
XX = Vector{Float64}[]
ZZ = Float64[]
PP = []
for i in eachindex(a)
    Z = quadgk(x -> exp(-ϕ4((x,), a[i])), -Inf, Inf)[1]
    push!(ZZ, Z)
    lf = My.LFProblem(1, ϕ4; dt = 0.05, nsteps = 100)
    @time X = first.(My.HMC(lf, a[i]))
    flush(stdout)
    push!(XX, X)
    P = plot()
    histogram!(X; norm=true, alpha=0.3, label="HMC LF sample", bin=100, c=i)
    plot!(x -> exp(-ϕ4(x, a[i]))/Z, -2, 2; label="exp(-ϕ2(x))/Z", lw=2, c=i)
    plot!(; legend=false, xtick=-2:0.5:2)
    title!("ϕ(x) = a(x² - 1)²  for  a = $(a[i])", titlefontsize=9)
    push!(PP, P)
end
plot(PP...; size=(800, 450), layout=(3, 2))

# %%
QQ = []
for i in eachindex(a)
    Q = plot(XX[i][1:10000]; ylim=(-1.5, 1.5), label="", c=i, lw=0.5)
    title!("ϕ(x) = a(x² - 1)²  for  a = $(a[i])", titlefontsize=9)
    push!(QQ, Q)
end
plot(QQ...; size=(800, 900), layout=(length(a), 1))

# %% [markdown]
# ## Baysian inference for a sample of the standard normal distribution

# %%
n = 10
sample_normal = randn(n)
f(y, m, s) = (y - m)^2/(2s^2) + log(s)
negloglik(w, sample) = sum(y -> f(y, w[1], exp(w[2])), sample)
lf = My.LFProblem(2, negloglik; dt = 0.1, nsteps = 30)

# %%
@time sample = My.HMC(lf, sample_normal; init = SVector(0.0, 0.0))
@time sample = My.HMC(lf, sample_normal; init = SVector(0.0, 0.0))
@time sample = My.HMC(lf, sample_normal; init = SVector(0.0, 0.0))

# %%
m, logs = first.(sample), last.(sample)
d = InterpKDE(kde((m, logs)))
x, y = range(extrema(m)...; length=201), range(extrema(logs)...; length=201)
heatmap(x, y, (x, y) -> pdf(d, x, y); size=(450, 400), xlabel="μ", ylabel="log(σ)")

# %% [markdown]
# ## Symbolics.jl example

# %%
dim = 2
@variables a[1:dim, 1:dim], x[1:dim]
aa, xx = collect.((a, x))

# %%
phi_sym = dot(xx, aa, xx)/2 |> expand |> simplify

# %%
dphi_sym = Symbolics.gradient(phi_sym, xx)

# %%
phi_rgf = build_function(phi_sym, xx, aa; expression = Val(false))
dphi_rgf = build_function(dphi_sym, xx, aa; expression = Val(false))[1]
lf = My.LFProblem(dim, phi_rgf, dphi_rgf);

# %%
param = A = @SMatrix [
     1  1/2
    1/2  1
]

# %%
@time sample = My.HMC(lf, param)
@time sample = My.HMC(lf, param)
@time sample = My.HMC(lf, param);

# %%
@btime My.HMC($lf, $param);

# %% [markdown]
# 自動微分を使うよりも計算が速くなっている.

# %%
X, Y = first.(sample), last.(sample)
d = InterpKDE(kde((X, Y)))
x, y = range(extrema(X)...; length=201), range(extrema(Y)...; length=201)
heatmap(x, y, (x, y) -> pdf(d, x, y); size=(450, 400), right_margin=3Plots.mm)

# %%
f(n) = mean(x -> x*x', @view sample[1:n])

n = 1:1000
S = f.(n)
S11 = (S -> S[1,1]).(S)
S22 = (S -> S[2,2]).(S)
S12 = (S -> S[1,2]).(S)

ymin = min(-1.5, minimum(S11), minimum(S22), minimum(S12))
ymax = max(2.5, maximum(S11), maximum(S22), maximum(S12))

plot(ylim = (ymin, ymax))
plot!(S11; label="s11", c=1)
hline!([inv(A)[1,1]]; label="", c=1, ls=:dash)
plot!(S22; label="s22", c=2)
hline!([inv(A)[2,2]]; label="", c=2, ls=:dash)
plot!(S12; label="s12", c=3)
hline!([inv(A)[1,2]]; label="", c=3, ls=:dash)

# %%
