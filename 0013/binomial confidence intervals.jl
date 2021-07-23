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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using RCall
macro sym_str(x) :(Symbol($(esc(x)))) end
@rimport binom as binom
@rimport exactci as exactci

# %%
using Distributions
using StatsPlots
using Roots
x ⪅ y = x < y || x ≈ y

pval_exact(dist, k) = sum(pdf(dist, j) for j in support(dist) if pdf(dist, j) ⪅ pdf(dist, k))
pval_exact(n, p, k) = pval_exact(Binomial(n, p), k)
pval_2x1sd(dist, k) = min(1, 2cdf(dist, k), 2ccdf(dist, k - 1))
pval_2x1sd(n, p, k) = pval_2x1sd(Binomial(n, p), k)
chisq_stat(n, p, k) = (k - n*p)^2/(n*p*(1 - p))
pval_chisq(n, p, k) = ccdf(Chisq(1), chisq_stat(n, p, k))
ci(pval, n, k, α=0.05) = find_zeros(p -> pval(n, p, k) - α, 0, 1)

# %%
n, k = 20, 3

# %%
binom.binom_confint(k, n)

# %%
@show exactci.binom_exact(k, n, tsmethod = "central")[sym"conf.int"][1:2]
@show exactci.binom_exact(k, n, tsmethod = "minlik")[sym"conf.int"][1:2]
@show exactci.binom_exact(k, n, tsmethod = "blaker")[sym"conf.int"][1:2];

# %%
@show ci(pval_exact, n, k)
@show ci(pval_2x1sd, n, k)
@show ci(pval_chisq, n, k);

# %%
n, k = 20, 3
plot(; xtick=0:0.1:1)
plot!(p -> pval_exact(n, p, k), 0, 1; label="exact two-sided (or minlik)")
plot!(p -> pval_2x1sd(n, p, k), 0, 1; label="exact doubled one-sided (or central)")
plot!(p -> pval_chisq(n, p, k), 0, 1; label="chi-squared (or Wilson)")

# %%
