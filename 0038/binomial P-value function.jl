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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# # 二項分布モデルでのP値函数のグラフ
#
# * 黒木玄
# * 2020-10-17
# $
# \newcommand\op{\operatorname}
# \newcommand\pvalue{\op{pvalue}}
# \newcommand\pdf{\op{pdf}}
# \newcommand\cdf{\op{cdf}}
# \newcommand\ccdf{\op{ccdf}}
# \newcommand\confint{\op{confint}}
# \newcommand\credint{\op{credint}}
# $

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#WilsonのP値函数の場合" data-toc-modified-id="WilsonのP値函数の場合-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>WilsonのP値函数の場合</a></span></li><li><span><a href="#WilsonのP値函数とClopper-PeasonのP値函数とベイズ版P値函数の比較" data-toc-modified-id="WilsonのP値函数とClopper-PeasonのP値函数とベイズ版P値函数の比較-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>WilsonのP値函数とClopper-PeasonのP値函数とベイズ版P値函数の比較</a></span></li><li><span><a href="#highest-density-interval-版のP値函数" data-toc-modified-id="highest-density-interval-版のP値函数-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>highest density interval 版のP値函数</a></span></li><li><span><a href="#信用区間にベイズ信用区間と同様の意味で「真の値」が含まれる確率" data-toc-modified-id="信用区間にベイズ信用区間と同様の意味で「真の値」が含まれる確率-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>信用区間にベイズ信用区間と同様の意味で「真の値」が含まれる確率</a></span></li></ul></div>

# %%
using DataFrames
using Distributions
using Roots
using StatsPlots
default(fmt=:png,
    titlefontsize=10, tickfontsize=6, guidefontsize=9, legendfontsize=9)

# %% [markdown]
# ## WilsonのP値函数の場合

# %%
function pvalue_bin_wilson(k, p; n=20)
    z = (k - n*p)/√(n*p*(1 - p))
    2ccdf(Normal(), abs(z))
end

# %%
n = 20
α = 0.05
x = range(-0.5, n+0.5, 500)
p = range(0, 1, 500)
heatmap(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n) ≥ α;
    clim=(0, 1), colorbar=false)
plot!(xguide="data k", yguide="parameter p")
plot!(xtick=0:n, ytick=0:0.1:1)
title!("Wilson's confidence intervals for n = $n, α = $α")
plot!(size=(580, 500))

# %%
n = 20
x = range(-0.5, n+0.5, 500)
p = range(0, 1, 500)
heatmap(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n); clim=(0, 1))
plot!(xguide="data k", yguide="parameter p")
plot!(xtick=0:n, ytick=0:0.1:1)
title!("Wilson's P-value function for n = $n")
plot!(size=(640, 500))

# %%
n, k = 20, 6
p = range(0, 1, 500)
plot(p, p -> pvalue_bin_wilson(k, p; n); label="")
plot!(xguide="parameter p", yguide="P-value")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
title!("Wilson's P-value function for n=$n, k=$k")

# %%
n = 20
x = range(-0.5, n+0.5, 211)
p = range(0, 1, 201)
@time anim = @animate for t in 0:3:359
    surface(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n); colorbar=false, camera=(t, 60))
    plot!(xguide="k", yguide="p", zguide="P-value")
    plot!(size=(640, 600), zlim=(0, 1))
end
gif(anim, "Wilson's P-value function n=$n (60).gif")

# %%
n = 20
x = range(-0.5, n+0.5, 211)
p = range(0, 1, 201)
@time anim = @animate for t in 0:3:359
    surface(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n); colorbar=false, camera=(t, 45))
    plot!(xguide="k", yguide="p", zguide="P-value")
    plot!(size=(640, 600), zlim=(0, 1))
end
gif(anim, "Wilson's P-value function n=$n (45).gif")

# %% [markdown]
# ## WilsonのP値函数とClopper-PeasonのP値函数とベイズ版P値函数の比較

# %%
safediv(x, y) = x==0 ? x : y==Inf ? zero(y) : x/y

# WilsonのP値函数
function pvalue_bin_wilson(k, p; n=20)
    z = safediv(k - n*p, √(n*p*(1 - p)))
    2ccdf(Normal(), abs(z))
end

# Clopper-PearsonのP値函数
function pvalue_bin_cp(k, p; n=20)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

# 事前分布Beta(γ, δ)によるBayes版のP値函数(equal-tailed interval版)
function pvalue_bin_bayesian_eti(k, p; n=20, γ=1, δ=1)
    beta = Beta(k+γ, n-k+δ)
    min(1, 2cdf(beta, p), 2ccdf(beta, p))
end

function plot_pvalue_functions(; n=20, k=6, γ=1, δ=1, kwargs...)
    p = range(0, 1, 500)
    plot(p, p -> pvalue_bin_wilson(k, p; n); label="Wilson")
    plot!(p, p -> pvalue_bin_cp(k, p; n); label="Clopper-Pearson", ls=:dash)
    plot!(p, p -> pvalue_bin_bayesian_eti(k, p; n, γ, δ); label="Bayesian (ETI)", ls=:dashdot)
    plot!(xguide="parameter p", yguide="P-value")
    plot!(xtick=0:0.1:1, ytick=0:0.1:1)
    title!("P-value functions for n=$n, k=$k, prior=Beta($γ, $δ)")
    plot!(; kwargs...)
end

# %%
plot_pvalue_functions(; n=20, k=6, γ=1, δ=1)

# %%
PP = []
for k in 0:11
    P = plot_pvalue_functions(; n=20, k, title="n=$n, k=$k", legend=false)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(4, 3))
plot!(legend=false)

# %%
plot_pvalue_functions(; n=20, k=6, γ=1//3, δ=1//3)

# %%
PP = []
for k in 0:11
    P = plot_pvalue_functions(; n=20, k, γ=1//3, δ=1//3, title="n=$n, k=$k", legend=false)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(4, 3))
plot!(legend=false)

# %% [markdown]
# ## highest density interval 版のP値函数
#
# 事後分布 $\op{Beta}(k+\gamma, n-k+\delta)$ の累積分布函数と確率密度函数をそれぞれ
#
# $$
# \Psi(p) = \cdf(\op{Beta}(k+\gamma, n-k+\delta), p), \quad
# \psi(p) = \pdf(\op{Beta}(k+\gamma, n-k+\delta), p)
# $$
#
# と書くとき, 仮説「成功確率は $p$ である」の highest density interval 版のP値が次のように定義される:
#
# $$
# \pvalue_{\op{Bayesian}\op{HDI}}(k|n,p,\gamma,\delta) = \Psi(p) + \Psi(p').
# $$
#
# ここで $p'$ は $\psi(p)=\psi(p')$ を満たす $p$ 以外のパラメータ値 $p'$ である.
#
# このP値の定義に対応する区間
#
# $$
# \credint_{\op{Bayesian}\op{HDI}}(k|n,\alpha,\gamma,\delta) =
# \{\, p\in[0,1] \mid \pvalue_{\op{Bayesian}\op{HDI}}(k|n,p,\gamma,\delta) \ge \alpha\,\}
# $$
#
# はhighest density interval 版の $100(1-\alpha)\%$ ベイズ信用区間, すなわち, 事後分布で測った確率が $1-\alpha$ になるようなパラメータの区間で長さが最短のものになる.

# %%
function pvalue_hdi(dist::ContinuousUnivariateDistribution, x₀; xlim = extrema(dist))
    p₀ = pdf(dist, x₀)
    m = mode(dist)
    f(x) = pdf(dist, x) - p₀
    if x₀ == m
        1.0
    elseif x₀ > m
        x₁ = find_zero(f, (xlim[begin], m))
        cdf(dist, x₁) + ccdf(dist, x₀)
    else
        x₁ = find_zero(f, (m, xlim[end]))
        cdf(dist, x₀) + ccdf(dist, x₁)
    end
end

# 事前分布Beta(γ, δ)によるBayes版のP値函数(highest density interval版)
function pvalue_bin_bayesian_hdi(k, p; n=20, γ=1, δ=1)
    beta = Beta(k+γ, n-k+δ)
    if k+γ ≤ 1
        return ccdf(beta, p)
    elseif n-k+δ ≤ 1
        return cdf(beta, p)
    end
    pvalue_hdi(beta, p)
end

function plot_pvalue_functions_hdi(; n=20, k=6, γ=1, δ=1, kwargs...)
    p = range(0, 1, 500)
    plot(p, p -> pvalue_bin_wilson(k, p; n); label="Wilson")
    plot!(p, p -> pvalue_bin_cp(k, p; n); label="Clopper-Pearson", ls=:dash)
    plot!(p, p -> pvalue_bin_bayesian_hdi(k, p; n, γ, δ); label="Bayesian (HDI)", ls=:dashdot)
    plot!(xguide="parameter p", yguide="P-value")
    plot!(xtick=0:0.1:1, ytick=0:0.1:1)
    title!("P-value functions for n=$n, k=$k, prior=Beta($γ, $δ)")
    plot!(; kwargs...)
end

# %%
plot_pvalue_functions_hdi(; n=20, k=6, γ=1, δ=1)

# %%
PP = []
for k in 0:11
    P = plot_pvalue_functions_hdi(; n=20, k, γ=1, δ=1, title="n=$n, k=$k", legend=false)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(4, 3))
plot!(legend=false)

# %% [markdown]
# ## 信用区間にベイズ信用区間と同様の意味で「真の値」が含まれる確率
#
# 結構頻繁に次のようなことを言う人をみかける:
#
# * 95%信頼区間に真の値が含まれる確率が95%になると考えるのは誤解であるが, ベイズ統計における95%信用区間であれば真の値が95%の確率で含まれると考えてよい.
#
# 事後分布にしたがってランダムに生成されたパラメータ値 $p$ が95%ベイズ信用区間に含まれる確率は95%になる.  これが「ベイズ統計における95%信用区間であれば真の値が95%の確率で含まれる」の正確な内容である.
#
# すなわち, 「ベイズ統計における95%信用区間であれば真の値が95%の確率で含まれる」における「真の値」は統計分析において真に興味がある現実における何かの値ではなく, 数学的フィクションである統計モデルの条件付き確率分布として定義される事後分布で生成された値に過ぎない.
#
# 「ベイズ統計における95%信用区間であれば真の値が95%の確率で含まれる」という発言を目にして, 「これこそ, 私が欲しかった性質である. 信頼区間よりもベイズ信用区間の方が解釈がわかりやすいし, 優れている!」のように思った人は, 過剰広告に騙されていることになる.
#
# 一様事前分布に関する事後分布については「事後分布にしたがってランダムに生成されたパラメータ値 $p$ が含まれる確率は95%になる」という性質は95%信頼区間についても近似的に成立している.  (注意: 実際には一様事前分布の場合に限らず, おとなしめの任意の事前分布についても同様である.)
#
# そのことを以下で数値的に確認しよう.

# %%
# Wilsonの信頼区間
function confint_bin_wilson(k; n=20, α=0.05)
    z = cquantile(Normal(0,1), α/2)
    p̂ = k/n
    a, b, c = 1 + z^2/n, p̂ + z^2/(2n), p̂^2
    sqrtD = √(b^2 - a*c)
    (b - sqrtD)/a, (b + sqrtD)/a
end

# 事後分布で測ったWilsonの信頼区間にpが含まれる確率
function prob_of_p_in_ci_wrt_posterior(k; n=20, α=0.05, γ=1, δ=1)
    beta = Beta(k+γ, n-k+δ)
    p_L, p_U = confint_bin_wilson(k; n, α)
    cdf(beta, p_U) - cdf(beta, p_L)
end

# 事後分布で測ったWilsonの信頼区間にpが含まれる確率のプロット
function plot_prob_of_p_in_ci_wrt_posterior(; n=20, γ=1, δ=1, kwargs...)
    α = 0.05
    ks = 0:n
    probs = prob_of_p_in_ci_wrt_posterior.(ks; n, α, γ, δ)

    plot(ks, probs; label="")
    hline!([0.95]; label="", ls=:dash)
    plot!(ylim=(0.9, 1.002))
    plot!(xtick=0:n, ytick=0:0.01:1)
    plot!(xguide="k", yguide="probability")
    title!("prob. of p ∈ 95% conf. int. w.r.t. posterior for n=$n, prior=Beta($γ, $δ)")
    plot!(; kwargs...)
end

# %%
plot_prob_of_p_in_ci_wrt_posterior(; n=20)

# %%
plot_prob_of_p_in_ci_wrt_posterior(; n=100, xtick=0:5:100)

# %%
plot_prob_of_p_in_ci_wrt_posterior(; n=1000, xtick=0:50:1000)

# %% [markdown]
# このように一様事前分布から得られる事後分布で測った95%信頼区間(Wilson)に含まれる確率は95%でよく近似されている.
#
# ただし, $k/n$ または $1-k/n$ が小さ過ぎる場合にはその確率は95%よりも大きめになることには注意が必要である.

# %%
