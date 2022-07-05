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
using RCall
using StatsPlots

default(fmt=:png, titlefontsize=10)
@rlibrary stats

function pvalue_clopper_pearson(model::DiscreteUnivariateDistribution, data)
    min(1, 2cdf(model, data), 2ccdf(model, data-1))
end

x ⪅ y = x < y || x ≈ y

function pvalue_sterne(model::DiscreteUnivariateDistribution, data)
    p0 = pdf(model, data)
    m = mode(model)
    pdf(model, m) ≈ p0 && return 1.0
    if data > m
        i = m - 1
        while !(pdf(model, i) ⪅ p0) i -= 1 end
        return cdf(model, i) + ccdf(model, data-1)
    else # data < m
        i = m + 1
        while !(pdf(model, i) ⪅ p0) i += 1 end
        return cdf(model, data) + ccdf(model, i-1)
    end
end

label = ["Sterne" "Clopper-Pearson" "R"]
ls = [:solid :dash :dashdot]
ytick = 0:0.05:1
xguide = "data"
yguide = "P-value"

# %%
n, p = 20, 0.5
model = Binomial(n, p)
datas = support(model)
pvals = [
    pvalue_sterne.(model, datas);;
    pvalue_clopper_pearson.(model, datas);;
    [rcopy(binom_test(x, n; p))[:p_value] for x in datas]
]
@show model datas
plot(datas, pvals; label, ls, xtick=datas, ytick, xguide, yguide)
title!("bionmial test for n=$n, p=$p")

# %%
n, p = 20, 0.37
model = Binomial(n, p)
datas = support(model)
pvals = [
    pvalue_sterne.(model, datas);;
    pvalue_clopper_pearson.(model, datas);;
    [rcopy(binom_test(x, n; p))[:p_value] for x in datas]
]
@show model datas
plot(datas, pvals; label, ls, xtick=datas, ytick, xguide, yguide)
title!("bionmial test for n=$n, p=$p")

# %%
m, n, r = 19, 21, 15
model = Hypergeometric(m, n, r)
datas = support(model)
matdatas = [[a m-a; r-a n-r+a] for a in datas]
pvals = [
    pvalue_sterne.(model, datas);;
    pvalue_clopper_pearson.(model, datas);;
    [rcopy(fisher_test(x))[:p_value] for x in matdatas]
]
@show model datas
plot(datas, pvals; label, ls, xtick=datas, ytick, xguide="data a", yguide)
title!("fisher test for a+b=$m, c+d=$n, a+c=$r")

# %%
λ = 5
model = Poisson(λ)
datas = 0:15
pvals = [
    pvalue_sterne.(model, datas);;
    pvalue_clopper_pearson.(model, datas);;
    [rcopy(poisson_test(x, λ))[:p_value] for x in datas]
]
@show model datas
plot(datas, pvals; label, ls, xtick=datas, ytick, xguide="data", yguide)
title!("Poisson test for λ=$λ")

# %%
