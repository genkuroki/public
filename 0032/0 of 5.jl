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
default(fmt=:png, titlefontsize=10)
safediv(x, y) = x==0 ? x : isinf(y) ? zero(y) : x/y
using Roots

# %%
function pvalue_wilson(n, k, p)
    μ = n*p
    sehat = √(n*p*(1-p))
    z = safediv(k - n*p, sehat)
    2ccdf(Normal(), abs(z))
end

# %%
pvalue_wilson(5, 0, 0.1)

# %%
find_zero(p -> pvalue_wilson(5, 0, p) - 0.05, 0.5)

# %%
function pvalue_clopper_pearson(n, k, p)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

# %%
pvalue_clopper_pearson(5, 0, 0.2)

# %%
find_zero(p -> pvalue_clopper_pearson(5, 0, p) - 0.05, 0.5)

# %%
function pvalue_bayes(n, k, p; a=0.5, b=0.5)
    beta = Beta(a+k, b+n-k)
     min(1, 2cdf(beta, p), 2ccdf(beta, p))
end

# %%
pvalue_bayes(5, 0, 0.1)

# %%
find_zero(p -> pvalue_bayes(5, 0, p) - 0.05, 0.5)

# %%
find_zero(p -> pvalue_bayes(5, 0, p; a=1, b=1) - 0.05, 0.5)

# %%
