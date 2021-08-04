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
#  http://www.snap-tck.com/room04/c01/stat/stat01/stat0106.html#note03

# %%
using Logging
disable_logging(Logging.Warn)

using Distributions
using Roots
using RCall
@rlibrary stats

x ⪅ y = x < y || x ≈ y

function pvalue(d, k, ::Val{:ts}) # :ts stands for "two-sides"
    sum(pdf(d, j) for j in support(d) if pdf(d, j) ⪅ pdf(d, k); init = 0.0)
end

function pvalue(d, k, ::Val{:dos}) # :dos stands for "doubled one-sided"
    min(1, 2cdf(d, k), 2ccdf(d, k-1))
end

function pvalue(a, b, c, d, ::Val{:fisher}; ω = 1.0)
    pvalue(FisherNoncentralHypergeometric(a+b, c+d, a+c, ω), a, Val(:ts))
end

function pvalue(a, b, c, d, ::Val{:fisher_dos}; ω = 1.0)
    pvalue(FisherNoncentralHypergeometric(a+b, c+d, a+c, ω), a, Val(:dos))
end

function confidence_interval(a, b, c, d, alg; α = 0.05)
    CI = exp.(find_zeros(t -> pvalue(a, b, c, d, alg; ω = exp(t)) - α, -20, 20))
    isone(length(CI)) ? CI[1] < 1 ? (0.0, CI[1]) : (CI[1], Inf) : (CI[1], CI[end])
end

function delta(a, b, c, d, ω=1.0)
    A = 1 - ω
    B = a + d + ω*(b + c)
    C = a*d - ω*b*c
    2*C/(-B - √(B^2 - 4*A*C))
end

_odds_ratio(a, b, c, d) = (ad = a*d; iszero(ad) ? ad : ad/(b*c))
function odds_ratio(a, b, c, d, ω=1.0)
    δ = delta(a, b, c, d, ω)
    _odds_ratio(a + δ, b - δ, c - δ, d + δ)
end

_chisq(a, b, c, d) = (a+b+c+d)*(a*d - b*c)^2/((a+b)*(c+d)*(a+c)*(b+d))
function chisq(a, b, c, d, ω=1.0)
    δ = delta(a, b, c, d, ω)
    iszero(δ) ? δ : δ^2 * (1/(a + δ) + 1/(b - δ) + 1/(c - δ) + 1/(d + δ))
end

function chisq_yates(a, b, c, d, ω=1.0)
    δ = delta(a, b, c, d, ω)
    m = max(0, abs(δ)-1/2)^2
    iszero(m) ? m : m * (1/(a + δ) + 1/(b - δ) + 1/(c - δ) + 1/(d + δ))
end

function pvalue(a, b, c, d, ::Val{:chisq}; ω = 1.0)
    ccdf(Chisq(1), chisq(a, b, c, d, ω))
end

function pvalue(a, b, c, d, ::Val{:chisq_yates}; ω = 1.0)
    ccdf(Chisq(1), chisq_yates(a, b, c, d, ω))
end

# %%
println("="^30 * " Test data:\n")

A = [16 4; 4 6]
@show A
println()

println("="^30 * " Calculated by R:\n")
@show rcopy(fisher_test(A))[:p_value]
@show rcopy(fisher_test(A))[:conf_int]
println()
@show rcopy(chisq_test(A, correct=false))[:p_value]
@show rcopy(chisq_test(A))[:p_value]

sleep(0.1)
println()

println("="^30 * " Calculated by Julia:\n")
@show pvalue(A..., Val(:fisher))
@show confidence_interval(A..., Val(:fisher))
@show pvalue(A..., Val(:fisher_dos))
@show confidence_interval(A..., Val(:fisher_dos))
println()
@show pvalue(A..., Val(:chisq))
@show confidence_interval(A..., Val(:chisq))
@show pvalue(A..., Val(:chisq_yates))
@show confidence_interval(A..., Val(:chisq_yates))
;

# %%
fisher_test(A)

# %%
chisq_test(A, correct = false)

# %%
chisq_test(A)

# %%
