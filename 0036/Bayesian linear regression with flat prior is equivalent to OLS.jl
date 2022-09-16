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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %%
ENV["COLUMNS"] = 120

using Distributions
using LinearAlgebra
normsquared(x) = dot(x, x)
using Random
Random.seed!(4649373)
using StatsPlots
using Turing

# %%
# Jeffreys事前分布などのimproper事前分布を定義するために以下が使われる.

"""
    PowerPos(p::Real)

The *positive power distribution* with real-valued parameter `p` is the improper distribution
of real numbers that has the improper probability density function

```math
f(x) = \\begin{cases}
0 & \\text{if } x \\leq 0, \\\\
x^p & \\text{otherwise}.
\\end{cases}
```
"""
struct PowerPos{T<:Real} <: ContinuousUnivariateDistribution
    p::T
end
PowerPos(p::Integer) = PowerPos(float(p))

Base.minimum(d::PowerPos{T}) where T = zero(T)
Base.maximum(d::PowerPos{T}) where T = T(Inf)

Base.rand(rng::Random.AbstractRNG, d::PowerPos) = rand(rng) + 0.5
function Distributions.logpdf(d::PowerPos, x::Real)
    T = float(eltype(x))
    return x ≤ 0 ? T(-Inf) : d.p*log(x)
end

Distributions.pdf(d::PowerPos, x::Real) = exp(logpdf(d, x))

# For vec support
function Distributions.loglikelihood(d::PowerPos, x::AbstractVector{<:Real})
    T = float(eltype(x))
    return any(xi ≤ 0 for xi in x) ? T(-Inf) : d.p*log(prod(x))
end

@doc PowerPos

# %%
@model function ols(x, y)
    σ² ~ PowerPos(-1)
    β₀ ~ Turing.Flat()
    β₁ ~ Turing.Flat()
    y ~ MvNormal(β₀ .+ β₁*x, σ²*I)
end

# %%
n = 20
t, e = 100randn(n), 100randn(n)
x = 10000 .+ t
y = 100000 .- 100t + 50e

X = x .^ (0:1)'
β̂₀, β̂₁ = X \ y
s² = 1/(n-2) * normsquared(y .- β̂₀ - β̂₁*x)
@show β̂₀ β̂₁ s²
scatter(x, y; label="data")

# %%
chn = sample(ols(x, y), NUTS(), MCMCThreads(), 2*10^4, 10);

# %%
chn

# %%
@show β̂₀ β̂₁ s²;

# %%
plot(chn; lw=0.5, alpha=0.5)

# %%
function posterior_marginal(x_new, x, y)
    n, r = length(x), 2
    X = x .^ (0:1)'
    β̂ = X \ y
    s² = normsquared(y - X*β̂)/(n-r)
    ρ = √(s²*(x_new'/(X'X)*x_new))
    x_new'β̂ + ρ * TDist(n-r)
end

function preddist_marginal(x_new, x, y)
    n, r = length(x), 2
    X = x .^ (0:1)'
    β̂ = X \ y
    s² = normsquared(y - X*β̂)/(n-r)
    ρ = √(s²*(1 + x_new'/(X'X)*x_new))
    x_new'β̂ + ρ * TDist(n-r)
end

# %%
post_β₀ = posterior_marginal([1, 0], x, y)
stephist(vec(chn[:β₀]); norm=true, label="MCMC")
plot!(post_β₀; label="theoretical", ls=:dash, lw=2)
title!("posterior of β₀")

# %%
post_β₁ = posterior_marginal([0, 1], x, y)
stephist(vec(chn[:β₁]); norm=true, label="MCMC")
plot!(post_β₁; label="theoretical", ls=:dash, lw=2)
title!("posterior of β₁")

# %%
