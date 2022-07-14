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
using Roots
using StatsFuns
using StatsPlots

default(fmt=:png, titlefontsize=10)
@rlibrary stats
safediv(x, y) = x==0 ? x : y==Inf ? zero(y) : x/y
x ⪅ y = x < y || x ≈ y

# %%
function pvalue_wilson(n, k, p)
    p̂ = k/n
    SE = √(p*(1-p)/n)
    z = safediv(p̂ - p, SE)
    2(1 - cdf(Normal(), abs(z)))
end

pvalue_wilson(30, 25, 2/3)

# %%
rcopy(prop_test(25, 30, p=2/3, correct=false))[:p_value]

# %%
function pvalue_wald(n, k, p)
    p̂ = k/n
    SEhat = √(p̂*(1-p̂)/n)
    z = safediv(p̂ - p, SEhat)
    2(1 - cdf(Normal(), abs(z)))
end

pvalue_wald(30, 25, 2/3)

# %%
function pvalue_logit_wald(n, k, p)
    p̂ = k/n
    SEhatinv = √(n*p̂*(1-p̂))
    z = (logit(p̂) - logit(p))*SEhatinv
    2(1 - cdf(Normal(), abs(z)))
end

pvalue_logit_wald(30, 25, 2/3)

# %%
function pvalue_bayes(n, k, p; conjprior=(0.5, 0.5))
    α, β = conjprior
    beta = Beta(α+k, β+n-k)
    min(1, 2cdf(beta, p), 2ccdf(beta, p))
end

pvalue_bayes(30, 25, 2/3)

# %%
function pvalue_clopper_pearson(n, k, p)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

pvalue_clopper_pearson(30, 25, 2/3)

# %%
function pvalue_sterne(n, k, p)
    bin = Binomial(n, p)
    p0 = pdf(bin, k)
    sum(pdf(bin, i) for i in support(bin) if pdf(bin, i) ⪅ p0)
end

pvalue_sterne(30, 25, 2/3)

# %%
rcopy(binom_test(25, 30, p=2/3))[:p_value]

# %%
function ci(pvaluefunc, n, k; α = 0.05)
    f(p) = pvaluefunc(n, k, p) - α
    find_zeros(f, 1e-4, 1 - 1e-4)
end

@show ci(pvalue_wilson, 30, 25)
@show ci(pvalue_wald, 30, 25)
@show ci(pvalue_logit_wald, 30, 25)
@show ci(pvalue_bayes, 30, 25)
@show ci(pvalue_clopper_pearson, 30, 25)
@show ci(pvalue_sterne, 30, 25)
@show rcopy(prop_test(25, 30, p=2/3, correct=false))[:conf_int]
@show rcopy(binom_test(25, 30, p=2/3))[:conf_int]
;

# %%
@show ci(pvalue_wilson, 30, 25; α=0.1)
@show ci(pvalue_wald, 30, 25; α=0.1)
@show ci(pvalue_logit_wald, 30, 25; α=0.1)
@show ci(pvalue_bayes, 30, 25; α=0.1)
@show ci(pvalue_clopper_pearson, 30, 25; α=0.1)
@show ci(pvalue_sterne, 30, 25; α=0.1)
@show rcopy(prop_test(25, 30, p=2/3, correct=false, var"conf.level"=0.9))[:conf_int]
@show rcopy(binom_test(25, 30, p=2/3, var"conf.level"=0.9))[:conf_int]
;

# %%
n, k = 30, 25

P1 = plot(p -> pvalue_wilson(n, k, p), 0.5, 1; label="Wilson")
plot!(p -> pvalue_wald(n, k, p); label="Wald", ls=:dash)
plot!(p -> pvalue_logit_wald(n, k, p); label="logit Wald", ls=:dashdot)
plot!(p -> pvalue_bayes(n, k, p); label="Bayes", ls=:dashdotdot)
plot!(legend=:topleft)
plot!(xguide="p", yguide="P-value")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
title!("data: n = $n, k = $k")

P2 = plot(p -> pvalue_clopper_pearson(n, k, p), 0.5, 1; label="Clopper-Pearson", c=4)
plot!(p -> pvalue_sterne(n, k, p); label="Sterne", ls=:dashdot, c=5)
plot!(legend=:topleft)
plot!(xguide="p", yguide="P-value")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
title!("data: n = $n, k = $k")

display(P1)
display(P2)

# %%
function power(pvaluefunc, n, p, q=p; α = 0.05)
    bin = Binomial(n, q)
    sum((pvaluefunc(n, k, p) < α)*pdf(bin, k) for k in support(bin))
end

@show power(pvalue_wilson, 30, 2/3, 25/30)
@show power(pvalue_wald, 30, 2/3, 25/30)
@show power(pvalue_logit_wald, 30, 2/3, 25/30)
@show power(pvalue_bayes, 30, 2/3, 25/30)
@show power(pvalue_clopper_pearson, 30, 2/3, 25/30)
@show power(pvalue_sterne, 30, 2/3, 25/30)
;

# %%
n, p = 30, 2/3

Q1 = plot(q -> power(pvalue_wilson, n, p, q), 0, 1; label="Wilson")
plot!(q -> power(pvalue_wald, n, p, q); label="Wald", ls=:dash)
plot!(q -> power(pvalue_logit_wald, n, p, q); label="logit Wald", ls=:dashdot)
plot!(q -> power(pvalue_bayes, n, p, q); label="Bayes", ls=:dashdotdot)
plot!(legend=:bottomleft)
plot!(xguide="alternative p", yguide="power")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
title!("n = $n, null p = $(round(p; digits=3))")

Q2 = plot(q -> power(pvalue_clopper_pearson, n, p, q), 0, 1; label="Clopper-Pearson", c=5)
plot!(q -> power(pvalue_sterne, n, p, q); label="Sterne", ls=:dash, c=6)
plot!(legend=:bottomleft)
plot!(xguide="alternative p", yguide="power")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
title!("n = $n, null p = $(round(p; digits=3))")

display(Q1)
display(Q2)

# %% [markdown]
# なるほど. logit Wald が使用されない理由がよくわかった!

# %%
n, p = 30, 2/3

Q1 = plot(q -> power(pvalue_wilson, n, p, q), 0, 1; label="Wilson")
plot!(q -> power(pvalue_wald, n, p, q); label="Wald", ls=:dash)
#plot!(q -> power(pvalue_logit_wald, n, p, q); label="logit Wald", ls=:dashdot)
plot!(q -> power(pvalue_bayes, n, p, q); label="Bayes", ls=:dashdotdot)
plot!(legend=:bottomleft)
plot!(xguide="alternative p", yguide="power")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
title!("n = $n, null p = $(round(p; digits=3))")

#Q2 = plot(q -> power(pvalue_clopper_pearson, n, p, q), 0, 1; label="Clopper-Pearson", c=5)
plot!(q -> power(pvalue_clopper_pearson, n, p, q), 0, 1; label="Clopper-Pearson", c=5)
plot!(q -> power(pvalue_sterne, n, p, q); label="Sterne", ls=:dash, c=6)
#plot!(legend=:bottomleft)
#plot!(xguide="alternative p", yguide="power")
#plot!(xtick=0:0.1:1, ytick=0:0.1:1)
#title!("n = $n, null p = $(round(p; digits=3))")

display(Q1)
#display(Q2)

# %% [markdown]
# WaldとBayes以外の3つの検出力は一致している.

# %%
n, p = 30, 2/3

Q1 = plot(q -> power(pvalue_wilson, n, p, q), 0, 1; label="Wilson")
#plot!(q -> power(pvalue_wald, n, p, q); label="Wald", ls=:dash)
#plot!(q -> power(pvalue_logit_wald, n, p, q); label="logit Wald", ls=:dashdot)
#plot!(q -> power(pvalue_bayes, n, p, q); label="Bayes", ls=:dashdotdot)
plot!(legend=:bottomleft)
plot!(xguide="alternative p", yguide="power")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
title!("n = $n, null p = $(round(p; digits=3))")

#Q2 = plot(q -> power(pvalue_clopper_pearson, n, p, q), 0, 1; label="Clopper-Pearson", c=5)
plot!(q -> power(pvalue_clopper_pearson, n, p, q), 0, 1; label="Clopper-Pearson", c=5)
plot!(q -> power(pvalue_sterne, n, p, q); label="Sterne", ls=:dash, c=6)
#plot!(legend=:bottomleft)
#plot!(xguide="alternative p", yguide="power")
#plot!(xtick=0:0.1:1, ytick=0:0.1:1)
#title!("n = $n, null p = $(round(p; digits=3))")

display(Q1)
#display(Q2)

# %%
binom_test(25, 30, p=19/30)

# %%
19/30

# %%
