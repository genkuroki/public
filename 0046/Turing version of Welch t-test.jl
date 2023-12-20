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
#     display_name: Julia 1.9.4
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using HypothesisTests

function student_t_test(X, Y; Δμ = 0.0)
    m, X̄, SX2 = length(X), mean(X), var(X)
    n, Ȳ, SY2 = length(Y), mean(Y), var(Y)
    S2 = ((m-1)*SX2 + (n-1)*SY2) / (m+n-2)
    sehat2 = S2 * (1/m + 1/n)
    tvalue = (X̄ - Ȳ - Δμ) / √sehat2
    df = m + n - 2
    pvalue = 2ccdf(TDist(df), abs(tvalue))
    (; pvalue, tvalue, sehat2, df)
end

pvalue_student_t_test(X, Y; Δμ = 0.0) = student_t_test(X, Y; Δμ).pvalue

function welch_t_test(X, Y; Δμ = 0.0)
    m, X̄, SX2 = length(X), mean(X), var(X)
    n, Ȳ, SY2 = length(Y), mean(Y), var(Y)
    sehat2 = SX2/m + SY2/n
    tvalue = (X̄ - Ȳ - Δμ) / √sehat2
    df = sehat2^2 / ((SX2/m)^2/(m-1) + (SY2/n)^2/(n-1))
    pvalue = 2ccdf(TDist(df), abs(tvalue))
    (; pvalue, tvalue, sehat2, df)
end

pvalue_welch_t_test(X, Y; Δμ = 0.0) =welch_t_test(X, Y; Δμ).pvalue

function confint_welch_t_test(X, Y; α = 0.05)
    m, X̄, SX2 = length(X), mean(X), var(X)
    n, Ȳ, SY2 = length(Y), mean(Y), var(Y)
    sehat2 = SX2/m + SY2/n
    df = sehat2^2 / ((SX2/m)^2/(m-1) + (SY2/n)^2/(n-1))
    c = cquantile(TDist(df), α/2)
    ci = [X̄ - Ȳ - c*√sehat2, X̄ - Ȳ + c*√sehat2]
    ci
end

function pvalue_mann_whitney_utest(X, Y)
    o = MannWhitneyUTest(X, Y)
    pvalue(o)
end

function brunner_munzel_test(X, Y; p=1/2)
    m, n = length(X), length(Y)
    phat = mean((x < y) + (x == y)/2 for x in X, y in Y)
    Hbarx = n*(1 - phat)
    Hbary = m*phat
    sx2 = 1/n^2 * 1/(m-1) * sum(x -> (sum((y < x) + (y == x)/2 for y in Y) - Hbarx)^2, X)
    sy2 = 1/m^2 * 1/(n-1) * sum(y -> (sum((x < y) + (x == y)/2 for x in X) - Hbary)^2, Y)
    sehat = √(sx2/m + sy2/n)
    tvalue = (phat - p)/sehat
    df = (sx2/m + sy2/n)^2 / ((sx2/m)^2/(m-1) + (sy2/n)^2/(n-1))
    pvalue = 2ccdf(TDist(df), abs(tvalue))
    (; phat, sehat, tvalue, df, pvalue, p)
end

pvalue_brunner_munzel_test(X, Y; p=1/2) = brunner_munzel_test(X, Y; p).pvalue

# %% [markdown]
# https://x.com/__fusion/status/1736851891068498229
#
# * 治療薬: 6.7, 6.8, 9.0, 8.8, 6.5, 7.7, 8.1, 7.5, 6.5, 9.5, 6.6, 5.9
# * 対照薬: 9.8, 9.6, 8.3, 8.9, 8.1, 7.9, 8.7, 8.3, 6.7, 8.8, 6.9, 10.0, 8.0, 7.8

# %%
Y = [6.7, 6.8, 9.0, 8.8, 6.5, 7.7, 8.1, 7.5, 6.5, 9.5, 6.6, 5.9]
X = [9.8, 9.6, 8.3, 8.9, 8.1, 7.9, 8.7, 8.3, 6.7, 8.8, 6.9, 10.0, 8.0, 7.8]

@show pvalue_student_t_test(X, Y)
@show pvalue_welch_t_test(X, Y)
@show pvalue_mann_whitney_utest(X, Y)
@show pvalue_brunner_munzel_test(X, Y)
@show confint_welch_t_test(X, Y);

# %%
using Distributions
using Random

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
ENV["COLUMNS"] = 200
using Turing
using LinearAlgebra
using StatsPlots
default(fmt=:png)

# %%
@model function model_FlatPos(X1, X0)
    sigma0 ~ FlatPos(0)
    sigma1 ~ FlatPos(0)
    mu0 ~ Flat()
    deltamu ~ Flat()
    mu1 = mu0 + deltamu
    X0 ~ MvNormal(fill(mu0, length(X0)), sigma0^2*I)
    X1 ~ MvNormal(fill(mu1, length(X1)), sigma1^2*I)
end

# %%
chn_FlatPos = sample(model_FlatPos(X, Y), NUTS(), MCMCThreads(), 10^5, 10);

# %%
chn_FlatPos

# %%
@show pvalue_welch_t_test(X, Y)
@show confint_welch_t_test(X, Y);

# %%
plot(chn_FlatPos[1:min(end, 10000)]; lw=0.5, alpha=0.5)

# %%
@model function model_PowerPos(X1, X0)
    sigma0 ~ PowerPos(-1)
    sigma1 ~ PowerPos(-1)
    mu0 ~ Flat()
    deltamu ~ Flat()
    mu1 = mu0 + deltamu
    X0 ~ MvNormal(fill(mu0, length(X0)), sigma0^2*I)
    X1 ~ MvNormal(fill(mu1, length(X1)), sigma1^2*I)
end

# %%
chn_PowerPos = sample(model_PowerPos(X, Y), NUTS(), MCMCThreads(), 10^5, 10);

# %%
chn_PowerPos

# %%
@show pvalue_welch_t_test(X, Y)
@show confint_welch_t_test(X, Y);

# %%
plot(chn_PowerPos[1:min(end, 10000)]; lw=0.5, alpha=0.5)

# %%
@model function model_Flat_logsigma(X1, X0)
    logsigma0 ~ Flat()
    logsigma1 ~ Flat()
    mu0 ~ Flat()
    deltamu ~ Flat()
    mu1 = mu0 + deltamu
    X0 ~ MvNormal(fill(mu0, length(X0)), exp(2logsigma0)*I)
    X1 ~ MvNormal(fill(mu1, length(X1)), exp(2logsigma1)*I)
end

# %%
chn_Flat_logsigma = sample(model_Flat_logsigma(X, Y), NUTS(), MCMCThreads(), 10^5, 10);

# %%
chn_Flat_logsigma

# %%
@show pvalue_welch_t_test(X, Y)
@show confint_welch_t_test(X, Y);

# %%
plot(chn_Flat_logsigma[1:min(end, 10000)]; lw=0.5, alpha=0.5)

# %%
