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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

# %%
function pvalue_ztest(x̄, n)
    2ccdf(Normal(), abs(√n * x̄))
end

function suffsize_ztest(μ; α=0.05, β=0.20)
    c = cquantile(Normal(), α/2)
    d = cquantile(Normal(), 1-β)
    n = ceil(Int, ((d - c)/μ)^2)
end

function power_ztest(μ, n; α=0.05)
    c = cquantile(Normal(), α/2)
    dist_x̄ = Normal(μ, 1/√n)
    cdf(dist_x̄, -c/√n) + ccdf(dist_x̄, c/√n)
end

# %%
@show μ = -0.5
@show n = suffsize_ztest(μ)
power_ztest(μ, n)

# %%
function sim(; α=0.05, β=0.02, n0=100, Niters=10^6)
    pval = zeros(Niters)
    for i in 1:Niters
        x̄0 = rand(Normal(0, 1/√n0))
        n = suffsize_ztest(x̄0; α, β)
        if n ≤ n0
            pval[i] = pvalue_ztest(x̄0, n0)
        else
            x̄1 = rand(Normal(0, 1/√(n-n0)))
            x̄ = (n0 * x̄0 + (n-n0) * x̄1) / n
            pval[i] = pvalue_ztest(x̄, n)
        end
    end
    @show mean(pval .< α)
    pval
end

# %%
sim(n0=100)

# %%
