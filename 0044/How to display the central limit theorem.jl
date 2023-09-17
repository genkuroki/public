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
using Distributions
using Random
using StatsBase: ecdf
using StatsPlots
default(fmt=:png)
using SymPy

function make_ecdf(A)
    _ecdf = ecdf(A)
    f(x) = _ecdf(x)
    f
end

ecdf_naive(A, x) = count(≤(x), A)/length(A)

distname(dist) = replace(string(dist), r"{[^\}]*}"=>"")

function distname(dist::MixtureModel)
    d = distname.(dist.components)
    d = replace.(d, r", σ=1.0"=>"")
    p = string.(dist.prior.p)
    name = "$(p[1]) $(d[1])"
    name *= prod(" + $(p[i]) $(d[i])" for i in 2:length(d))
    name
end

# %% [markdown]
# $\newcommand\op{\operatorname}$
#
# 期待値 $0$, 分散 $1$ の一様分布 $\op{Uniform}(-\sqrt{3}, \sqrt{3})$ のモーメント母函数は
#
# $$
# M(t) = \frac{1}{2\sqrt{3}}\int_{-\sqrt{3}}^{\sqrt{3}} e^{tx}\,dx
# = \frac{e^{\sqrt{3}\,t} - e^{-\sqrt{3}\,t}}{2\sqrt{3}\,t}
# $$
#
# キュムラント母函数 $\log M(t)$ を求めよう.

# %%
@vars t
expr = log((exp(√Sym(3)*t) - exp(-√Sym(3)*t))/(2√Sym(3)*t))
Eq(expr, expr.series(t; n=12).simplify())

# %% [markdown]
# 離散一様分布 $\op{DiscreteUniform}(1, N)$ の期待値と分散とモーメント母函数はそれぞれ
#
# $$
# \mu = \frac{N+1}{2}, \quad
# \sigma^2 = \frac{N^2-1}{12}, \quad
# M(t) = \frac{1}{N}\sum_{k=1}^N e^{tk}
# = \frac{e^{(N+1)t} - e^t}{N(e^t - 1)}
# = \frac{e^{Nt}-1}{N(1-e^{-t})}.
# $$
#
# 標準化されたキュムラント母函数
#
# $$
# \log M(t/\sigma) - \frac{\mu}{\sigma} t
# $$
#
# を $N=6$ の場合について求めよう.

# %%
@vars t
N = Sym(6)
μ = (N + 1)/2
σ = √((N^2 - 1)/12)
expr = (exp(N*t) - 1)/(N*(1 - exp(-t)))
expr = log(expr(t => t/σ)) - μ/σ*t
Eq(expr, expr.series(t; n=12).simplify())

# %%
function ecdf_samplemean(dist, n; L=10^6)
    X̄ = Vector{Float64}(undef, L)
    Xtmp = [Vector{eltype(dist)}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(dist, Xtmp[tid])
        X̄[i] = mean(X)
    end
    make_ecdf(X̄)
end

# %%
@show dist = Uniform()
@show n = 12
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
X̄ = f._ecdf.sorted_values
histogram(X̄; norm=true, alpha=0.3, bin=50, label="")
plot!(Normal(μ, σ/√n); label="", lw=1.5)

# %%
@show dist = Uniform()
@show n = 12
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
X̄ = f._ecdf.sorted_values
histogram(X̄; norm=true, alpha=0.3, bin=200, label="")
plot!(Normal(μ, σ/√n); label="", lw=1.5)

# %%
@show dist = DiscreteUniform(1, 6)
@show n = 12
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
X̄ = f._ecdf.sorted_values
histogram(X̄; norm=true, alpha=0.3, bin=1-1/24:1/12:6+1/24, label="")
plot!(Normal(μ, σ/√n); label="", lw=1.5)

# %%
@show dist = DiscreteUniform(1, 6)
@show n = 12
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
X̄ = f._ecdf.sorted_values
histogram(X̄; norm=true, alpha=0.3, bin=1-1/24:1/60:6+1/24, label="")
plot!(Normal(μ, σ/√n); label="", lw=1.5)

# %%
@show dist = MixtureModel([Normal(), Normal(40)], [1-1/40, 1/40])
@show n = 100
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
X̄ = f._ecdf.sorted_values
histogram(X̄; norm=true, alpha=0.3, bin=20, label="")
plot!(Normal(μ, σ/√n); label="", lw=1.5)

# %%
@show dist = MixtureModel([Normal(), Normal(40)], [1-1/40, 1/40])
@show n = 100
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
X̄ = f._ecdf.sorted_values
histogram(X̄; norm=true, alpha=0.3, bin=200, label="")
plot!(Normal(μ, σ/√n); label="", lw=1.5)

# %%
@show dist = Uniform()
@show n = 12
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
plot(z -> f(μ + σ*z/√n), -4, 4; label="ecdf of standardized sample mean")
plot!(z -> cdf(Normal(), z), -4, 4; ls=:dot, label="cdf of standardized normal")
plot!(xguide="z", yguide="cumulative probability")
title!("$(distname(dist)), n=$n", titlefontsize=10)

# %%
@show dist = DiscreteUniform(1, 6)
@show n = 12
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
plot(z -> f(μ + σ*z/√n), -4, 4; label="ecdf of standardized sample mean")
plot!(z -> cdf(Normal(), z), -4, 4; ls=:dot, label="cdf of standard normal")
plot!(xguide="z", yguide="cumulative probability")
title!("$(distname(dist)), n=$n", titlefontsize=10)

# %%
@show dist = MixtureModel([Normal(), Normal(40)], [1-1/40, 1/40])
@show n = 100
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
z = range(-4, 4, 1000)
plot(z, z -> f(μ + σ*z/√n); label="ecdf of standardized sample mean")
plot!(z, z -> cdf(Normal(), z); ls=:dot, label="cdf of standard normal")
plot!(xguide="z", yguide="cumulative probability")
title!("$(distname(dist)), n=$n", titlefontsize=10)

# %%
@show dist = Uniform()
@show n = 12
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
z = range(-4, 4, 1000)
plot(cdf.(Normal(), z), @.(f(μ + σ*z/√n)); label="")
plot!(identity, 0, 1; label="", ls=:dot)
plot!(xguide="cdf of standard normal", yguide="ecdf of standardized sample mean")
title!("$(distname(dist)), n=$n", titlefontsize=10)
plot!(size=(400, 400))

# %% tags=[]
@show dist = DiscreteUniform(1, 6)
@show n = 12
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
z = range(-4, 4, 1000)
plot(cdf.(Normal(), z), @.(f(μ + σ*z/√n)); label="")
plot!(identity, 0, 1; label="", ls=:dot)
plot!(xguide="cdf of standard normal", yguide="ecdf of standardized sample mean")
title!("$(distname(dist)), n=$n", titlefontsize=10)
plot!(size=(400, 400))

# %%
@show dist = MixtureModel([Normal(), Normal(40)], [1-1/40, 1/40])
@show n = 100
@show μ, σ = mean(dist), std(dist)
f = ecdf_samplemean(dist, n)
z = range(-4, 4, 1000)
plot(cdf.(Normal(), z), @.(f(μ + σ*z/√n)); label="")
plot!(identity, 0, 1; label="", ls=:dot)
plot!(xguide="cdf of standard normal", yguide="ecdf of standardized sample mean")
title!("$(distname(dist)), n=$n", titlefontsize=8)
plot!(size=(400, 400))

# %%
