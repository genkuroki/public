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
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

# %%
function pvalue_clopper_pearson(k, n, p)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

pvalue_clopper_pearson(6, 10, 1/4)

# %%
2cdf(Beta(6, 10-6+1), 1/4)

# %%
x ⪅ y = x < y || x ≈ y

function pvalue_sterne(k, n, p)
    bin = Binomial(n, p)
    sum(pdf(bin, j) for j in support(bin) if pdf(bin, j) ⪅ pdf(bin, k))
end

pvalue_sterne(6, 10, 1/4)

# %%
function pvalue_wilson(k, n, p)
    bin = Binomial(n, p)
    μ, σ = mean(bin), std(bin)
    normal = Normal(μ, σ)
    min(1, 2cdf(normal, k), 2ccdf(normal, k))
end

pvalue_wilson(6, 10, 1/4)

# %%
function pvalue_wald(k, n, p)
    bin = Binomial(n, p)
    p̂ = k/n
    μ, σ = mean(bin), √(n * p̂ * (1 - p̂))
    normal = Normal(μ, σ)
    min(1, 2cdf(normal, k), 2ccdf(normal, k))
end

pvalue_wald(6, 10, 1/4)

# %%
