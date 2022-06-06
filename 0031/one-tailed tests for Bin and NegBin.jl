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

# %% [markdown]
# # 二項分布と負の二項分布の片側検定のP値
#
# * 黒木玄
# * 2022-06-06
#
# $
# \newcommand\op{\operatorname}
# \newcommand\Beta{\op{Beta}}
# \newcommand\Binomial{\op{Binomial}}
# \newcommand\NegativeBinomial{\op{NegativeBinomial}}
# \newcommand\Bin{\op{Bin}}
# \newcommand\NegBin{\op{NegBin}}
# \newcommand\pmf{\op{pmf}}
# \newcommand\pdf{\op{pdf}}
# \newcommand\cdf{\op{cdf}}
# \newcommand\ccdf{\op{ccdf}}
# $

# %% [markdown]
# 一般に次が成立している:
#
# $$
# \begin{aligned}
# &
# \frac{\int_0^p t^{k-1}(1-t)^{n-k}\,dt}{B(k, n-k+1)} =
# \sum_{i=k}^n \binom{n}{i} p^i (1-p)^{n-i} =
# \sum_{j=k}^n \binom{j-1}{k-1} p^k (1-p)^{j-k},
# \\ &
# \frac{\int_p^1 t^{k}(1-t)^{n-k-1}\,dt}{B(k+1, n-k)} =
# \sum_{i=0}^k \binom{n}{i} p^i (1-p)^{n-i} =
# \sum_{j=n+1}^\infty \binom{j-1}{(k+1)-1} p^{k+1} (1-p)^{j-(k+1)},
# \\ &
# \frac{\int_p^1 t^{k-1}(1-t)^{n-k-1}\,dt}{B(k, n-k)} =
# \sum_{i=0}^{k-1} \binom{n-1}{i} p^i (1-p)^{(n-1)-i} =
# \sum_{j=n}^\infty \binom{j-1}{k-1} p^k (1-p)^{j-k}.
# \end{aligned}
# $$
#
# これは「試行回数 $n$, 成功回数 $k$」というデータについて, 以下のようになることを含む:
#
# $$
# \begin{aligned}
# &
# (\text{二項分布 $\Binomial(n, p)$ での仮説 $p\le p_0$ の片側検定のP値}) \\ &=
# (\text{負の二項分布 $\NegativeBinomial(k, p)$ での仮説 $p\le p_0$ の片側検定のP値}) \\ &=
# (\text{ベータ分布 $\Beta(k, n-k+1)$ 内で仮説 $p\le p_0$ が成立する確率})
# \\ &
# \quad
# \\ &
# (\text{二項分布 $\Binomial(n, p)$ での仮説 $p\ge p_0$ の片側検定のP値}) \\ &=
# (\text{ベータ分布 $\Beta(k+1, n-k)$ 内で仮説 $p\ge p_0$ が成立する確率})
# \\ &
# \quad
# \\ &
# (\text{負の二項分布 $\NegativeBinomial(k, p)$ での仮説 $p\ge p_0$ の片側検定のP値}) \\ &=
# (\text{ベータ分布 $\Beta(k, n-k)$ 内で仮説 $p\ge p_0$ が成立する確率}).
# \end{aligned}
# $$
#
# これは以下が成立していることを意味している:
#
# (1) 二項分布 $\Binomial(n, p)$ と負の二項分布 $\NegativeBinomial(k, p)$ のどちらにおいても, 仮説 $p\le p_0$ の片側検定のP値は, improper事前分布 $\Beta(0, 1)$ から定まる事後分布内でその仮説が成立する確率に等しい.
#
# (2) 二項分布 $\Binomial(n, p)$ での仮説 $p\ge p_0$ の片側検定のP値は, improper事前分布 $\Beta(1, 0)$ から定まる事後分布内でその仮説が成立する確率に等しい.
#
# (3) 負の二項分布 $\NegativeBinomial(k, p)$ での仮説 $p\ge p_0$ の片側検定のP値は, improper事前分布 $\Beta(0, 0)$ から定まる事後分布内で仮説 $p\ge p_0$ が成立する確率に等しい.
#
# このとき, (2)と(3)が一致していないことに注意せよ.

# %%
using Distributions
#using StatsPlots

# %%
n, k, p = 12, 3, 0.5
@show n k p
println()

@show sum(binomial(n, i)*p^i*(1-p)^(n-i) for i in k:n)
@show sum(binomial(j-1, k-1)*p^k*(1-p)^(j-k) for j in k:n)
@show cdf(Beta(k, n-k+1), p)
println()

@show sum(binomial(n, i)*p^i*(1-p)^(n-i) for i in 0:k)
@show 1 - sum(binomial(j-1, (k+1)-1)*p^(k+1)*(1-p)^(j-(k+1)) for j in k+1:n)
@show ccdf(Beta(k+1, n-k), p)
println()

@show sum(binomial(n-1, i)*p^i*(1-p)^((n-1)-i) for i in 0:k-1)
@show 1 - sum(binomial(j-1, k-1)*p^k*(1-p)^(j-k) for j in k:n-1)
@show ccdf(Beta(k, n-k), p)
;

# %%
n, k, p = 12, 3, 0.5
@show n k p
println()

@show ccdf(Binomial(n, p), k-1)
@show cdf(NegativeBinomial(k, p), n-k)
@show cdf(Beta(k, n-k+1), p)
println()

@show cdf(Binomial(n, p), k)
@show ccdf(NegativeBinomial(k+1, p), (n+1)-1-(k+1))
@show ccdf(Beta(k+1, n-k), p)
println()

@show cdf(Binomial(n-1, p), k-1)
@show ccdf(NegativeBinomial(k, p), n-1-k)
@show ccdf(Beta(k, n-k), p)
;

# %%
n, k, p = 24, 7, 0.5
@show n k p
println()

@show ccdf(Binomial(n, p), k-1)
@show cdf(NegativeBinomial(k, p), n-k)
@show cdf(Beta(k, n-k+1), p)
println()

@show cdf(Binomial(n, p), k)
@show ccdf(NegativeBinomial(k+1, p), (n+1)-1-(k+1))
@show ccdf(Beta(k+1, n-k), p)
println()

@show cdf(Binomial(n-1, p), k-1)
@show ccdf(NegativeBinomial(k, p), n-1-k)
@show ccdf(Beta(k, n-k), p)
;

# %%
