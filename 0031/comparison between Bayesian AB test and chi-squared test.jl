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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Distributions
using StatsPlots
using StatsBase

# %%
using RCall

@rimport base as R
R.Sys_setenv(LANG = "en")

@rimport stats as stats
@rlibrary exact2x2
@rlibrary exactci

# %%
using SymPy: @vars, expand, sympy
@vars a b c d ω δ positive=true
expr = (a-δ)*(d-δ) - ω*(b+δ)*(c+δ)
sympy.Poly(expand(expr), δ)

# %% tags=[]
function pval_chisq(a, b, c, d, ω)
    A = 1 - ω
    B = a + d + ω*(b + c)
    C = a*d - ω*b*c
    δ = 2C/(B + √(B^2 - 4A*C))
    χ² = δ^2 * (1/(a-δ) + 1/(b+δ) + 1/(c+δ) + 1/(d-δ))
    ccdf(Chisq(1), χ²)
end

# %%
α, β = 0.5, 0.5
a, b, c, d = 5, 20, 11, 14
beta1 = Beta(α + a, β + b)
beta2 = Beta(α + c, β + d)

L = 10^6
R1 = rand(beta1, L)
R2 = rand(beta2, L)

OR = @. R1/(1-R1)/(R2/(1-R2))
ecdf_OR = ecdf(OR)
pval_bayes(x) = min(1, 2ecdf_OR(x), 2 - 2ecdf_OR(x))

P1 = stephist(OR; norm=true, xlim=(0.0, 1.5), label="posteriot of OR")
P2 = plot(x -> pval_bayes(x), 0.0, 1.5; label="Bayesian P-value")
plot!(x -> pval_chisq(a, b, c, d, x); label="χ² P-value", ls=:dash)

@show [a b; c d]
@show pval_bayes(1.0)
@show pval_chisq(a, b, c, d, 1.0)

plot(P1, P2; size=(800, 250))

# %%
exact2x2([a b; c d]; plot=true)

# %%
exact2x2([a b; c d]; plot=true, tsmethod="central")

# %%
stats.chisq_test([a b; c d]; correct=true)

# %%
stats.chisq_test([a b; c d]; correct=false)

# %%
