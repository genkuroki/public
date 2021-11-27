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
#     display_name: Julia 1.6.4
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# # 通常の信頼区間とベイズ版信用区間の比較

# %%
using Distributions
using Roots
using Plots
using StatsFuns: logistic, logit

# \approx TAB → ≈
# \lessapprox TAB → ⪅
x ⪅ y = x < y || x ≈ y

# %%
# Bernoulli分布モデルでの通常の信頼区間

"""P値函数"""
function pvalue(n, p, k)
    bin = Binomial(n, p)
    p0 = pdf(bin, k)
    s = sum(pdf(bin, j) for j in support(bin) if pdf(bin, j) ⪅ p0)
    min(1, s)
end

"""信頼区間函数"""
function confidence_interval(n, k; α = 0.05)
    f(t) = pvalue(n, logistic(t), k) - α
    CI = logistic.(find_zeros(f, -10.0, 10.0))
    if length(CI) < 2
        return 2k ≤ n ? [0.0, first(CI)] : [first(CI), 1.0]
    else
        return [first(CI), last(CI)]
    end
end

# %%
# Bernoulli分布モデルと共役事前分布に関するベイズ版信用区間

"""共役事後分布"""
posterior_dist(n, k; a = 0.5, b = 0.5) = Beta(a + k, b + n - k)

"""ベイズ版P値函数"""
function pvalue_bayes(n, p, k; a = 0.5, b = 0.5)
    posterior = posterior_dist(n, k; a, b)
    v0 = logpdf(posterior, p)
    f(t) = logpdf(posterior, logistic(t)) - v0
    m = params(posterior) |> ((α, β),) -> (α - 1)/(α + β - 2)
    if m ≤ 0
        s = ccdf(posterior, p)
    elseif m ≥ 1
        s = cdf(posterior, p)
    elseif p ≤ m
        q = logistic(find_zero(f, min(50, logit(m) + 1)))
        s = cdf(posterior, p) + ccdf(posterior, q)
    else
        q = logistic(find_zero(f, max(-50, logit(m) - 1)))
        s = cdf(posterior, q) + ccdf(posterior, p)
    end
    min(1, s)
end

"""ベイズ版信用区間函数"""
function credible_interval(n, k; α = 0.05, a = 0.5, b = 0.5)
    g(t) = pvalue_bayes(n, logistic(t), k; a, b) - α
    posterior = posterior_dist(n, k; a, b)
    m = params(posterior) |> ((α, β),) -> (α - 1)/(α + β - 2)
    if m ≤ 0
        p = logistic.(find_zero(g, 0.0))
        return [0.0, p]
    elseif m ≥ 1
        p = logistic.(find_zero(g, 0.0))
        return [p, 1.0]
    else
        L, R = max(-50, logit(m) - 5), min(50, logit(m) + 5)
        CI = logistic.(find_zeros(g, L, R))
        if length(CI) < 2
            return 2k ≤ n ? [0.0, first(CI)] : [first(CI), 1.0]
        else
            return [first(CI), last(CI)]
        end
    end
end

# %%
# 信頼区間とベイズ版信用区間の比較

function plot_cis(n, k; α = 0.05, a = 0.5, b = 0.5)
    confint = confidence_interval(n, k; α)
    credint = credible_interval(n, k; α, a, b)

    bin = Binomial(n, k/n)
    μ, σ = mean(bin)/n, std(bin)/n
    p = range(max(1e-8, μ - 3.5σ), min(1-1e-8, μ + 4σ); length=400)
    plot(p, p -> pvalue_bayes(n, p, k); label="Bayesian p-value function", c = 1, lw=1.5)
    plot!(p, p -> pvalue(n, p, k); label="ordinary p-value function", c = 2, ls=:dash, lw=1.5)
    plot!(; ytick=0:0.05:1)
    plot!(credint, fill(1.1α, 2); label="Bayesian credible interval", c=:blue, lw=2)
    plot!(confint, fill(0.9α, 2); label="confidence interval", c=:red, ls=:dash, lw=2)
    title!("n = $n, k = $k, α = $α, prior = Beta($a, $b)"; titlefontsize=10)
end

# %%
plot_cis(10, 3)

# %%
plot_cis(100, 30)

# %%
plot_cis(1000, 300)

# %% [markdown]
# $n$ を大きくするにつれて通常の信頼区間とベイズ版信用区間は数値的によく一致するようになる。
#
# このようなことはシンプルなモデルでは一般的に広く成立している.

# %%
