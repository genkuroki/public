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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using Roots

function pvalue_score(k, n, p)
    phat = k/n
    z = (phat - p) / √(p*(1-p)/n)
    2ccdf(Normal(), abs(z))
end

function confint_score(k, n, α=0.05)
    f(p) = pvalue_score(k, n, p) - α
    find_zeros(f, (eps(), 1-eps()))
end

# %%
n, p = 10000, 1/6
k = 0.170n
@show pvalue_score(k, n, p)
@show confint_score(k, n);

# %%
n, p = 50000, 1/6
k = 0.170n
@show pvalue_score(k, n, p)
@show confint_score(k, n);

# %%
function expecval(f, n, p; m = 5)
    bin = Binomial(n, p)
    μ, σ = mean(bin), std(bin)
    kmin = max(0, round(Int, μ-m*σ))
    kmax = min(n, round(Int, μ+m*σ))
    sum(k -> f(k) * pdf(bin, k), kmin:kmax)
end

# %%
n, p₀, p₁, α = 24_500, 1/6, 0.160, 0.05
expecval(k -> pvalue_score(k, n, p₀) < α, n, p₁)

# %%
n, p₀, p₁, α = 100_000, 1/6, 0.170, 0.05
expecval(k -> pvalue_score(k, n, p₀) < α, n, p₁)

# %%
n, p₀, p₁, α = 400_000, 1/6, 0.165, 0.05
expecval(k -> pvalue_score(k, n, p₀) < α, n, p₁)

# %%
n, p₀, p₁, α = 620_000, 1/6, 0.168, 0.05
expecval(k -> pvalue_score(k, n, p₀) < α, n, p₁)

# %%
n, p₀, p₁, α = 2_500_000, 1/6, 0.166, 0.05
expecval(k -> pvalue_score(k, n, p₀) < α, n, p₁)

# %%
