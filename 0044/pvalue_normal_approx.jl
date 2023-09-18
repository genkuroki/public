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
using StatsPlots
default(fmt=:png)

function pvalue_normal_approx(nulldist, x)
    m, s = mean(nulldist), std(nulldist)
    z = (x - m)/s
    2ccdf(Normal(), abs(z))
end

make_ecdf(Y) = (_ecdf = ecdf(Y); f(x) = _ecdf(x))

function plot_ecdf_pval(F_pval)
    plot(F_pval, 0, 1; label="")
    plot!(identity, 0, 1; label="", ls=:dot, c=:black, alpha=0.5)
    plot!(xtick=0:0.1:1, ytick=0:0.1:1)
    plot!(xguide="nominal significance level α",
        yguide="probability of P-value ≤ α")
    plot!(size=(400, 400))
end

@show pvalue_normal_approx(Binomial(100, 0.5), 60);

@show nulldist = Binomial(100, 1/3)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
@show nulldist = Binomial(100, 1/3)
@show altdist = Binomial(100, 1/2)
X = rand(altdist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
plot(p -> pvalue_normal_approx(Binomial(100, p), 30), 0, 1; label="data: n=100, x=30")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="parameter p", yguide="P-value")

# %%
@show nulldist = Poisson(30)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
@show nulldist = NegativeBinomial(30, 0.7)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
@show nulldist = Hypergeometric(200, 200, 200)
X = rand(nulldist, 10^6)
pval = pvalue_normal_approx.(nulldist, X)
F_pval = make_ecdf(pval)

println("probability of P-value ≤ 5% = ", F_pval(0.05))

plot_ecdf_pval(F_pval)

# %%
