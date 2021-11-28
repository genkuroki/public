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

"""P値函数(exact版)"""
function pvalue_exact(n, p, k)
    bin = Binomial(n, p)
    p0 = pdf(bin, k)
    s = sum(pdf(bin, j) for j in support(bin) if pdf(bin, j) ⪅ p0)
    min(1, s)
end

"""P値函数(正規分布近似)"""
function pvalue_normal(n, p, k)
    bin = Binomial(n, p)
    normal = Normal(mean(bin), std(bin))
    min(1, 2cdf(normal, k), 2ccdf(normal, k))
end

"""信頼区間函数"""
function confidence_interval(pvalue_func, n, k; α = 0.05)
    f(t) = pvalue_func(n, logistic(t), k) - α
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
posterior_dist(n, k; a = 1, b = a) = Beta(a + k, b + n - k)

"""ベイズ版P値函数"""
function pvalue_bayes(n, p, k; a = 1, b = a)
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
function credible_interval(n, k; α = 0.05, a = 1, b = a)
    confidence_interval((n, p, k) -> pvalue_bayes(n, p, k; a, b), n, k; α)
end

# %%
# 信頼区間とベイズ版信用区間の比較

function plot_cis(n, k; α = 0.05, a = 1, b = a)
    ci_exact = confidence_interval(pvalue_exact, n, k; α)
    ci_normal = confidence_interval(pvalue_normal, n, k; α)
    ci_bayes = credible_interval(n, k; α, a, b)

    bin = Binomial(n, k/n)
    μ, σ = mean(bin)/n, std(bin)/n
    p = range(max(1e-8, μ - 3.5σ), min(1-1e-8, μ + 4σ); length=400)
    plot(p, p -> pvalue_bayes(n, p, k; a, b); label="Bayesian p-value", c = 1, lw=1.5)
    plot!(p, p -> pvalue_normal(n, p, k); label="p-value (normal approx.)", c = 2, ls=:dash, lw=1.5)
    plot!(p, p -> pvalue_exact(n, p, k); label="p-value (exact)", c = 3, ls=:dashdot, lw=1.5)
    plot!(; ytick=0:0.05:1)
    plot!(ci_bayes, fill(1.2α, 2); label="Bayesian credible interval", c=1, lw=2)
    plot!(ci_normal, fill(1.0α, 2); label="conf. int. (normal approx.)", c=2, ls=:dash, lw=2)
    plot!(ci_exact, fill(0.8α, 2); label="confidence interval (exact)", c=3, ls=:dashdot, lw=2)
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
plot_cis(1000, 30)

# %%
plot_cis(1000, 3)

# %%
plot_cis(1000, 3; a = 0.5)

# %%
plot_cis(1000, 3; a = 0)

# %%
