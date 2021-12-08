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
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %%
using Distributions
using StatsPlots

# %% [markdown]
# 負の二項分布の確率函数は
#
# $$
# P(k) = \binom{k+r-1}{k} p^r (1 - p)^k \quad (k = 0,1,2,\ldots)
# $$
#
# であり, $k$ は成功確率 $p$ のBernoulli試行で初めて $r$ 回成功するまでの失敗する回数を意味している. 
#
# これのパラメータを $r = \alpha$, $p = 1/(N\theta)$ とおくと, $N$ が大きいとき $k/N$ は近似的に分布 $\operatorname{Gamma}(\alpha, \theta)$ に従う.

# %%
function plot_scaled_negbin(; α = 3, θ = 5, N = 100)
    gamma = Gamma(α, θ)
    negbin = LocationScale(0, 1/N, NegativeBinomial(α, 1/(N*θ)))
    
    xmax = quantile(gamma, 0.995)
    plot(; legend=:bottomright)
    plot!(x -> cdf(negbin, x), 0, xmax; label="NegativeBinomial(r = α, p = 1/(Nθ)) scaled by 1/N")
    plot!(x -> cdf(gamma, x), 0, xmax; label="Gamma(α = $α, θ = $θ)", ls=:dash)
    title!("cdf of scaled NegativeBinomial distribution for N = $N"; titlefontsize=10)
end

# %%
plot_scaled_negbin(N = 1)

# %%
plot_scaled_negbin(N = 10)

# %%
plot_scaled_negbin(N = 100)

# %%
α = 3
θ = 5
N = 100
gamma = Gamma(α, θ)
negbin = LocationScale(0, 1/N, NegativeBinomial(α, 1/(N*θ)))
X = rand(negbin, 10^5)

xmax = 60
stephist(X; norm=true, bin=0:xmax, label="NegativeBinomial(r = α, p = 1/(Nθ)) scaled by 1/N")
plot!(gamma, 0, xmax; label="Gamma(α = $α, θ = $θ)", ls=:dash)
title!("sample of scaled NegativeBinomial distribution for N = $N"; titlefontsize=10)

# %% [markdown]
# 二項分布の確率函数は
#
# $$
# P(k) = \binom{n}{k} p^k (1 - p)^{n-k} \quad (k = 0,1,\ldots,n)
# $$
#
# であり, $k$ は成功確率 $p$ の $n$ 回のBernoulli試行中の成功回数を意味している.
#
# これのパラメータを $p = \lambda/n$ とおくと, $n$ が大きいとき $k$ は近似的に分布 $\operatorname{Poisson}(\lambda)$ に従う.

# %%
function plot_bin(; λ = 5, n = 100)
    poisson = Poisson(λ)
    bin(n) = Binomial(n, λ/n)
    
    xmax = quantile(poisson, 0.995)
    plot(; legend=:bottomright)
    plot!(x -> cdf(bin(n), x), 0, xmax; label="Binomial(n, p = λ/n)")
    plot!(x -> cdf(poisson, x), 0, xmax; label="Poisson(λ = $λ)", ls=:dash)
    title!("cdf of Binomial distribution for n = $n"; titlefontsize=10)
end

# %%
plot_bin(n = 5)

# %%
plot_bin(n = 10)

# %%
plot_bin(n = 100)

# %%
