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
using StatsBase: ecdf
make_ecdf(Y) = (_ecdf = ecdf(Y); f(x) = _ecdf(x))
using StatsPlots
default(fmt=:png)

function pvalue_normal_approx(nulldist, x)
    m, s = mean(nulldist), std(nulldist)
    z = (x - m)/s
    2ccdf(Normal(), abs(z))
end

@show pvalue_normal_approx(Binomial(100, 0.5), 60);

@show nulldist = Binomial(100, 1/3)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F(0.05))

plot([F identity], 0, 1; label="", ls=[:solid :dash], size=(400, 400))

# %%
@show nulldist = Binomial(100, 1/3)
@show altdist = Binomial(100, 1/2)
X = rand(altdist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F(0.05))

plot([F identity], 0, 1; label="", ls=[:solid :dash], size=(400, 400))

# %%
@show nulldist = Poisson(30)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F(0.05))

plot([F identity], 0, 1; label="", ls=[:solid :dash], size=(400, 400))

# %%
@show nulldist = NegativeBinomial(30, 0.7)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F(0.05))

plot([F identity], 0, 1; label="", ls=[:solid :dash], size=(400, 400))

# %%
@show nulldist = Hypergeometric(200, 200, 200)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F(0.05))

plot([F identity], 0, 1; label="", ls=[:solid :dash], size=(400, 400))

# %%
