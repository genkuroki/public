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
default(fmt=:png, titlefontsize=10)

safediv(x, y) = x == 0 ? zero(x/y) : x/y
x ⪅ y = x < y || x ≈ y

function pvalue_score(k, n, p)
    phat = k/n
    sehat = √(p * (1 - p) / n)
    z = safediv(phat - p, sehat)
    2ccdf(Normal(), abs(z))
end

function pvalue_clopper_pearson(k, n, p)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k - 1))
end

function pvalue_sterne(k, n, p)
    bin = Binomial(n, p)
    sum(pdf(bin, j) for j in support(bin) if pdf(bin, j) ⪅ pdf(bin, k))
end

function expectval(f, dist; m=6)
    kmin, kmax = extrema(dist)
    μ, σ = mean(dist), std(dist)
    kmin = max(kmin, round(Int, μ - m*σ))
    kmax = min(kmax, round(Int, μ + m*σ))
    sum(k -> f(k) * pdf(dist, k), kmin:kmax)
end

function power(pvaluefunc, n, p0, p1; α=0.05)
    expectval(k -> pvaluefunc(k, n, p0) < α, Binomial(n, p1))
end

# %%
p0, p1 = 0.4, 0.6
ns = 30:70
plot()
plot!(ns, n -> power(pvalue_score, n, p0, p1); label="pvalue_score", c=1)
plot!(ns, n -> power(pvalue_clopper_pearson, n, p0, p1); label="pvalue_clopper_pearson", ls=:dashdot, c=2)
plot!(ns, n -> power(pvalue_sterne, n, p0, p1); label="pvalue_sterne", ls=:dash, c=3)
plot!(xguide="n", yguide="power", ytick=0:0.1:1, legend=:bottomright)
title!("power of binomial tests for null=$p0 and alt=$p1")

# %%
p0, p1 = 0.4, 0.6
ns = 30:70
plot()
plot!(ns, n -> power(pvalue_score, n, p0, p1); label="pvalue_score", c=1)
plot!(ns, n -> power(pvalue_clopper_pearson, n, p0, p1); label="pvalue_clopper_pearson", ls=:dashdot, c=2)
#plot!(ns, n -> power(pvalue_sterne, n, p0, p1); label="pvalue_sterne", ls=:dash, c=3)
plot!(xguide="n", yguide="power", ytick=0:0.1:1, legend=:bottomright)
title!("power of binomial tests for null=$p0 and alt=$p1")

# %%
p0, p1 = 0.4, 0.6
ns = 30:70
plot()
#plot!(ns, n -> power(pvalue_score, n, p0, p1); label="pvalue_score", c=1)
plot!(ns, n -> power(pvalue_clopper_pearson, n, p0, p1); label="pvalue_clopper_pearson", ls=:dashdot, c=2)
plot!(ns, n -> power(pvalue_sterne, n, p0, p1); label="pvalue_sterne", ls=:dash, c=3)
plot!(xguide="n", yguide="power", ytick=0:0.1:1, legend=:bottomright)
title!("power of binomial tests for null=$p0 and alt=$p1")

# %%
p0, p1 = 0.4, 0.6
ns = 30:70
plot()
plot!(ns, n -> power(pvalue_score, n, p0, p1); label="pvalue_score", c=1)
#plot!(ns, n -> power(pvalue_clopper_pearson, n, p0, p1); label="pvalue_clopper_pearson", ls=:dashdot, c=2)
plot!(ns, n -> power(pvalue_sterne, n, p0, p1); label="pvalue_sterne", ls=:dash, c=3)
plot!(xguide="n", yguide="power", ytick=0:0.1:1, legend=:bottomright)
title!("power of binomial tests for null=$p0 and alt=$p1")

# %%
p0, p1 = 0.1, 0.2
ns = 40:140
plot()
plot!(ns, n -> power(pvalue_score, n, p0, p1); label="pvalue_score", c=1)
plot!(ns, n -> power(pvalue_clopper_pearson, n, p0, p1); label="pvalue_clopper_pearson", ls=:dashdot, c=2)
plot!(ns, n -> power(pvalue_sterne, n, p0, p1); label="pvalue_sterne", ls=:dash, c=3)
plot!(xguide="n", yguide="power", ytick=0:0.1:1, legend=:bottomright)
title!("power of binomial tests for null=$p0 and alt=$p1")

# %%
p0, p1 = 0.1, 0.2
ns = 40:140
plot()
plot!(ns, n -> power(pvalue_score, n, p0, p1); label="pvalue_score", c=1)
plot!(ns, n -> power(pvalue_clopper_pearson, n, p0, p1); label="pvalue_clopper_pearson", ls=:dashdot, c=2)
#plot!(ns, n -> power(pvalue_sterne, n, p0, p1); label="pvalue_sterne", ls=:dash, c=3)
plot!(xguide="n", yguide="power", ytick=0:0.1:1, legend=:bottomright)
title!("power of binomial tests for null=$p0 and alt=$p1")

# %%
p0, p1 = 0.1, 0.2
ns = 40:140
plot()
#plot!(ns, n -> power(pvalue_score, n, p0, p1); label="pvalue_score", c=1)
plot!(ns, n -> power(pvalue_clopper_pearson, n, p0, p1); label="pvalue_clopper_pearson", ls=:dashdot, c=2)
plot!(ns, n -> power(pvalue_sterne, n, p0, p1); label="pvalue_sterne", ls=:dash, c=3)
plot!(xguide="n", yguide="power", ytick=0:0.1:1, legend=:bottomright)
title!("power of binomial tests for null=$p0 and alt=$p1")

# %%
p0, p1 = 0.1, 0.2
ns = 40:140
plot()
plot!(ns, n -> power(pvalue_score, n, p0, p1); label="pvalue_score", c=1)
#plot!(ns, n -> power(pvalue_clopper_pearson, n, p0, p1); label="pvalue_clopper_pearson", ls=:dashdot, c=2)
plot!(ns, n -> power(pvalue_sterne, n, p0, p1); label="pvalue_sterne", ls=:dash, c=3)
plot!(xguide="n", yguide="power", ytick=0:0.1:1, legend=:bottomright)
title!("power of binomial tests for null=$p0 and alt=$p1")

# %%
p0, p1 = 0.1, 0.15
ns = 280:400
plot()
plot!(ns, n -> power(pvalue_score, n, p0, p1); label="pvalue_score", c=1)
plot!(ns, n -> power(pvalue_clopper_pearson, n, p0, p1); label="pvalue_clopper_pearson", ls=:dashdot, c=2)
plot!(ns, n -> power(pvalue_sterne, n, p0, p1); label="pvalue_sterne", ls=:dash, c=3)
plot!(xguide="n", yguide="power", ytick=0:0.05:1, legend=:bottomright)
title!("power of binomial tests for null=$p0 and alt=$p1")

# %%
