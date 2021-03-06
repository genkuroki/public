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

# %% [markdown]
# # Clopper-Pearsonの信頼区間
#
# * 黒木玄
# * 2021-12-05
#
# $\newcommand\on{\operatorname}$
#
# 関連ノートブック
#
# * [止め方で結果が変わる？](https://nbviewer.org/github/genkuroki/public/blob/main/0025/%E6%AD%A2%E3%82%81%E6%96%B9%E3%81%A7%E7%B5%90%E6%9E%9C%E3%81%8C%E5%A4%89%E3%82%8F%E3%82%8B%EF%BC%9F.ipynb)
# * [ベイズハッキング](https://nbviewer.org/github/genkuroki/public/blob/main/0025/Bayes%20hacking.ipynb)

# %%
using Distributions
using StatsPlots
using Memoization
using Roots
using StatsFuns: logistic, logit

# %% [markdown]
# ## Clopper-Pearsonの信頼区間のP値函数の様々な表示
#
# ### 二項分布とベータ分布の関係
#
# 右辺の部分積分の繰り返しによって以下の公式を示すことができる(詳しい計算については「[止め方で結果が変わる？](https://nbviewer.org/github/genkuroki/public/blob/main/0025/%E6%AD%A2%E3%82%81%E6%96%B9%E3%81%A7%E7%B5%90%E6%9E%9C%E3%81%8C%E5%A4%89%E3%82%8F%E3%82%8B%EF%BC%9F.ipynb)」を参照):
#
# $$
# \begin{aligned}
# &
# \sum_{j=0}^k \binom{n}{j} \theta^j (1 - \theta)^{n-j} =
# \frac
# {\int_\theta^1 t^{(k+1)-1} (1 - t)^{(n-k)-1}\,dt}
# {B(k+1, n-k)} =
# \on{ccdf}(\on{Beta}(k+1, n-k), \theta),
# \\ &
# \sum_{j=k}^n \binom{n}{j} \theta^j (1 - \theta)^{n-j} =
# \frac
# {\int_0^\theta t^{k-1} (1 - t)^{(n-k+1)-1}\,dt}
# {B(k, n-k+1)} =
# \on{cdf}(\on{Beta}(k, n-k+1), \theta).
# \end{aligned}
# $$
#
# 前者の公式の左辺は二項分布 $\on{Binomial}(n, \theta)$ において成功回数が $k$ 以下になる確率であり, 「パラメータの値は $\theta$ 以上である」という帰無仮説の片側検定の $P$ 値とみなされる. 前者の公式の右辺はベイズ統計での対応する場合にimproper共役事前分布 $\on{Beta}(1, 0)$ から得られる事後分布で測ったパラメータの値が $\theta$ 以上になる確率になっている. 後者の公式は前者から導かれる. 
#
# これらの公式は片側検定に関する通常のP値とimproper共役事前分布を用いたベイズ統計の事後分布で測った帰無仮説が成立する確率がぴったり等しいことを意味している.

# %% [markdown]
# ### ベータ分布を使ったClopper-Pearsonの信頼区間のP値函数の表示
#
# 上の結果を使うと二項分布モデルでの片側検定のP値の2倍で定義した両側検定のP値
#
# $$
# \on{pvalue\_dos}(n, k, \theta) =
# \min\left\{
# 1,\;
# 2\sum_{j=0}^k \binom{n}{j} \theta^j (1 - \theta)^{n-j},\;
# 2\sum_{j=k}^n \binom{n}{j} \theta^j (1 - \theta)^{n-j}
# \right\}
# $$
#
# をベータ分布の $\on{cdf}$, $\on{ccdf}$ で表すことができる:
#
# $$
# \begin{aligned}
# \on{pvalue\_dos}(n, k, \theta) &=
# \min\left\{
# 1,\;
# 2\frac
# {\int_\theta^1 t^{(k+1)-1} (1 - t)^{(n-k)-1}\,dt}
# {B(k+1, n-k)},\;
# 2\frac
# {\int_0^\theta t^{k-1} (1 - t)^{(n-k+1)-1}\,dt}
# {B(k, n-k+1)}
# \right\}
# \\ &=
# \min\{
# 1,\;
# 2\on{ccdf}(\on{Beta}(k+1, n-k), \theta),\;
# 2\on{cdf}(\on{Beta}(k, n-k+1), \theta)
# \}.
# \end{aligned}
# $$
#
# ここで $\on{dos}$ は doubuled one-side の略である. この $\on{pvalue\_dos}(n, k, \theta)$ は __Clopper-Pearsonの信頼区間__ に対応するP値函数になっている.

# %% [markdown]
# ### F分布とベータ分布の関係
#
# 自由度 $\nu_1, \nu_2$ のF分布の確率密度函数の $dx$ 倍は以下のように書ける:
#
# $$
# \on{pdf}(\on{FDist}(\nu_1, \nu_2), x)\,dx = 
# \frac{1}{B(\nu_1/2, \nu_2/2)}
# \left(\frac{\nu_1 x}{\nu_1 x + \nu_2}\right)^{\nu_1/2}
# \left(1 - \frac{\nu_1 x}{\nu_1 x + \nu_2}\right)^{\nu_2/2}
# x^{-1}
# \,dx \quad (x > 0).
# $$
#
# これを
#
# $$
# t = \frac{\nu_1 x}{\nu_1 x + \nu_2}, \quad
# x = \frac{\nu_2}{\nu_1}\frac{t}{1 - t}, \quad 0 < t < 1
# $$
#
# で変数変換すると,
#
# $$
# \begin{aligned}
# \on{pdf}(\on{FDist}(\nu_1, \nu_2), x)\,dx &=
# \frac{1}{B(\nu_1/2, \nu_2/2)}
# t^{\nu_1/2}
# (1 - t)^{\nu_2/2}
# \,\frac{\nu_1}{\nu_2}\frac{1-t}{t}
# \,\frac{\nu_2}{\nu_1}\frac{dt}{(1-t)^2}
# \\ &=
# \frac{1}{B(\nu_1/2, \nu_2/2)} t^{\nu_1/2 - 1} (1 - t)^{\nu_2/2 - 1} \,dt
# \\ &=
# \on{pdf}(\on{Beta}(\nu_1/2, \nu_2/2), t)\,dt.
# \end{aligned}
# $$
#
# これは上の変数変換によって, 自由度 $\nu_1, \nu_2$ のF分布がパラメータ $\nu_1/2, \nu_2/2$ のベータ分布に変換されることを意味している. ゆえに
#
# $$
# \theta = \frac{\alpha x}{\alpha x + \beta}, \quad
# x = \frac{\beta}{\alpha}\frac{\theta}{1 - \theta}
# $$
#
# とおくと,
#
# $$
# \on{cdf}(\on{Beta}(\alpha, \beta), \theta) = \on{cdf}(\on{FDist}(2\alpha, 2\beta), x), \quad
# \on{ccdf}(\on{Beta}(\alpha, \beta), \theta) = \on{ccdf}(\on{FDist}(2\alpha, 2\beta), x).
# $$

# %% [markdown]
# ### F分布を使ったClopper-Pearsonの信頼区間のP値函数の表示
#
# 以上の結果を使えば __Clopper-Pearsonの信頼区間__ に対応するP値函数をF分布を使って計算することもできる:
#
# $$
# \on{pvalue\_dos}(n, k, \theta) =
# \min\begin{Bmatrix}
# 1 \\
# 2\on{ccdf}\left(\on{FDist}(2(k+1), 2(n-k)), \frac{n-k}{k+1}\frac{\theta}{1-\theta}\right) \\
# 2\on{cdf}\left(\on{FDist}(2k, 2(n-k+1)), \frac{n-k+1}{k}\frac{\theta}{1-\theta}\right)
# \end{Bmatrix}.
# $$

# %%
@memoize function pvalue_dos_naive(n, k, θ)
    bin = Binomial(n, θ)
    min(1, 2sum(pdf(bin, j) for j in 0:k), 2sum(pdf(bin, j) for j in k:n))
end

@memoize function pvalue_dos(n, k, θ)
    k ≤ 0 && return min(1, 2ccdf(Beta(k+1, n-k), θ))
    k ≥ n && return min(1, 2cdf(Beta(k, n-k+1), θ))
    min(1, 2ccdf(Beta(k+1, n-k), θ), 2cdf(Beta(k, n-k+1), θ))
end

@memoize function pvalue_dos_fdist(n, k, θ)
    k ≤ 0 && return min(1, 2ccdf(FDist(2(k+1), 2(n-k)), (n-k)/(k+1)*θ/(1-θ)))
    k ≥ n && return min(1, 2cdf(FDist(2k, 2(n-k+1)), (n-k+1)/k*θ/(1-θ)))
    min(1,
        2ccdf(FDist(2(k+1), 2(n-k)), (n-k)/(k+1)*θ/(1-θ)),
        2cdf(FDist(2k, 2(n-k+1)), (n-k+1)/k*θ/(1-θ)))
end

# %%
# Clopper-Pearsonの信頼区間に対応する3つのP値函数が一致することの確認

n = 10
k = 0:n
θ = 0:0.1:1

p1 = pvalue_dos_naive.(n, k, θ')

# %%
p2 = pvalue_dos.(n, k, θ')

# %%
p1 .≈ p2

# %%
p3 = pvalue_dos_fdist.(n, k, θ')

# %%
p3 .≈ p2

# %% [markdown]
# ### 負の二項分布とベータ分布の関係
#
# 上と同様に右辺の部分積分の繰り返しによって次の公式も示すことができる:
#
# $$
# \sum_{m=n}^\infty \binom{m-1}{k-1} \theta^k (1 - \theta)^{m-k} =
# \frac
# {\int_\theta^1 t^{k-1} (1 - t)^{(n-k)-1}\, dt}
# {B(k, n-k)}
# \quad (n \ge k \ge 1,\ 0 < \theta \le 1).
# $$
#
# 左辺の $\binom{m-1}{k-1} \theta^k (1 - \theta)^{m-k}$ は成功確率を意味するパラメータの値が $\theta$ のBernoulli試行をちょうど $k$ 回成功するまで繰り返したときの試行回数がちょうど $m$ 回になる確率なので, 左辺の和はちょうど $k$ 回の成功するまでの試行回数が $n$ 回以上になる確率になっており, 「パラメータの値は $\theta$ 以上である」という帰無仮説の片側検定の $P$ 値とみなされる.  右辺はベイズ統計での対応する場合にimproper事前分布 $\on{Beta}(0, 0)$ から得られる事後分布で測ったパラメータの値が $\theta$ 以上になる確率になっている.
#
# この公式も片側検定に関する通常のP値とimproper共役事前分布を用いたベイズ統計の事後分布で測った帰無仮説が成立する確率がぴったり等しくなっていることを意味している.

# %% [markdown]
# `Distributions.jl` における負の二項分布の定義は以下の通り:
#
# $$
# \on{pdf}(\on{NegativeBinomial}(r, \theta), k) = \binom{k+r-1}{k} p^r (1 - p)^k \quad (k=0,1,2,\ldots).
# $$
#
# ゆえに
#
# $$
# \on{pdf}(\on{NegativeBinomial}(k, \theta), m-k) =
# \binom{m-1}{m-k} \theta^k (1 - \theta)^{m-k} =
# \binom{m-1}{k-1} \theta^k (1 - \theta)^{m-k} \quad (m\ge k).
# $$

# %%
@memoize pvalue_negbin_naive(n, k, θ) =
    1 - sum(pdf(NegativeBinomial(k, θ), m - k) for m in k:n-1)

@memoize pvalue_negbin(n, k, θ) = ccdf(Beta(k, n-k), θ)

# %%
# 負の二項分布モデルでの片側検定に関する2つのP値函数が一致することの確認

k = 5
n = k+1:15
θ = 0.1:0.1:1

p1 = pvalue_negbin_naive.(n, k, θ')

# %%
p2 = pvalue_negbin.(n, k, θ')

# %%
p1 .≈ p2

# %% [markdown]
# ### もう１つのP値函数
#
# 上の片側確率の2倍版と異なる二項分布の両側検定のP値を次のように定めることもできる:
#
# $$
# \on{pvalue\_exact}(n, k, \theta) =
# \sum_{0\le j\le n,\; P(n,j,\theta)\le P(n,k,\theta)} P(n, k, \theta).
# $$
#
# ここで $
# P(n, k, \theta) = \binom{n}{j} \theta^j (1 - \theta)^{n-j}
# $.

# %%
x ⪅ y = x < y || x ≈ y # 以下の計算で ≤使うと失敗する

@memoize function pvalue_exact(n, k, θ)
    bin = Binomial(n, θ)
    p0 = pdf(bin, k)
    sum(pdf(bin, j) for j in support(bin) if pdf(bin, j) ⪅ p0)
end

# %%
# 別のP値函数が定義されたことの確認

n = 4
k = 0:n
θ = 0:0.1:1

p1 = pvalue_exact.(n, k', θ)
p2 = pvalue_dos.(n, k', θ)
p1 - p2

# %% [markdown]
# ## 信頼区間
#
# ### P値函数を用いた一般的な信頼区間の定義の仕方
#
# これらのＰ値函数 $\on{pvalue\_func} = \on{pvalue\_dos}, \on{pvalue\_exact}$ から信頼度 $1-\alpha$ の信用区間を次のように定義できる:
#
# $$
# \on{confidence\_interval}(\on{pvalue\_func}, n, k, \alpha) =
# \{\, \theta \in [0, 1] \mid \on{pvalue\_func}(n, k, \theta) \ge \alpha \,\}.
# $$
#
# これは $\on{pvalue\_func}(n, k, \theta) = \alpha$ となる $\theta$ を得ることによって数値計算できる.

# %%
@memoize function confidence_interval(pvalue_func, n, k, α = 0.05)
    f(t) = pvalue_func(n, k, logistic(t)) - α
    ci = logistic.(find_zeros(f, -10, 10))
    length(ci) < 2 && return 2k ≤ n ? (0.0, first(ci)) : (first(ci), 1.0) 
    (first(ci), last(ci))
end

# %%
[(k, confidence_interval(pvalue_dos, 10, k)) for k in 0:10]

# %%
[(k, confidence_interval(pvalue_exact, 10, k)) for k in 0:10]

# %% [markdown]
# ### Clopper-Pearsonの信頼区間のベータ分布やF分布を用いた表示
#
# 特に $\on{pvalue\_func} = \on{pvalue\_dos}$ の場合には上に述べたことから, 以下が成立することがわかる:
#
# $$
# \on{confidence\_interval}(\on{pvalue\_dos}, n, k, \alpha) = [\theta_L, \theta_U]
# $$
#
# ここで
#
# $$
# \begin{aligned}
# &
# \on{ccdf}(\on{Beta}(k+1, n-k), \theta_U) =
# \frac
# {\int_{\theta_U}^1 t^{(k+1)-1} (1 - t)^{(n-k)-1}\,dt}
# {B(k+1, n-k)} = \frac{\alpha}{2},
# \\ &
# \on{cdf}(\on{Beta}(k, n-k+1), \theta_L) =
# \frac
# {\int_0^{\theta_L} t^{k-1} (1 - t)^{(n-k+1)-1}\,dt}
# {B(k, n-k+1)} = \frac{\alpha}{2}.
# \end{aligned}
# $$
#
# これらを満たす $\theta_L, \theta_U$ はBeta分布に関する $\on{cquantile}$, $\on{quantile}$ 函数で計算できる:
#
# $$
# \theta_U = \on{cquantile}(\on{Beta}(k+1, n-k), \alpha/2), \quad
# \theta_L = \on{quantile}(\on{Beta}(k, n-k+1), \alpha/2).
# $$
#
# $\on{confidence\_interval}(\on{pvalue\_dos}, n, k, \alpha)$ を __Clopper-Pearsonの信頼区間__ と呼ぶ.

# %% [markdown]
# Clopper-Pearsonの信頼区間をF分布を用いて表示することもできる:
#
# $$
# \begin{aligned}
# &
# \theta_U = \frac{(k+1)x_U}{(k+1)x_U + (n-k)}, \quad
# x_U = \on{cquantile}(\on{FDist}(2(k+1), 2(n-k)), \alpha/2),
# \\ &
# \theta_L = \frac{k x_L}{k x_L + (n-k+1)}, \quad
# x_L = \on{quantile}(\on{FDist}(2k, 2(n-k+1)), \alpha/2).
# \end{aligned}
# $$

# %%
@memoize function confidence_interval_dos(n, k, α = 0.05)
    θ_U = k ≥ n ? 1.0 : cquantile(Beta(k+1, n-k), α/2)
    θ_L = k ≤ 0 ? 0.0 :  quantile(Beta(k, n-k+1), α/2)
    (θ_L, θ_U)
end

@memoize function confidence_interval_dos_fdist(n, k, α = 0.05)
    x_U = k ≥ n ? 1.0 : cquantile(FDist(2(k+1), 2(n-k)), α/2)
    x_L = k ≤ 0 ? 0.0 :  quantile(FDist(2k, 2(n-k+1)), α/2)
    θ_U = (k+1)*x_U/((k+1)*x_U + (n-k))
    θ_L = k*x_L/(k*x_L + (n-k+1))
    (θ_L, θ_U)
end

# %%
# confidence_interval(pvalue_dos, n, k) と confidence_interval_dos(n, k) が等しいことの確認

[(k, confidence_interval(pvalue_dos, 10, k) .≈ confidence_interval_dos(10, k)) for k in 0:10]

# %%
# confidence_interval(pvalue_dos, n, k) と confidence_interval_dos_fdist(n, k) が等しいことの確認

[(k, confidence_interval(pvalue_dos, 10, k) .≈ confidence_interval_dos_fdist(10, k)) for k in 0:10]

# %%
# 多重ディスパッチで効率化 (`find_zeros` 函数による計算は効率が悪い)

confidence_interval(::typeof(pvalue_dos), n, k, α = 0.05) = confidence_interval_dos(n, k, α)

# %%
methods(confidence_interval)

# %%
# 多重ディスパッチによる効率化の確認

@which confidence_interval(pvalue_dos, 10, 1)

# %% [markdown]
# ## P値函数と信頼区間のプロット
#
# ### 与えられたデータに対するP値函数と信頼区間のプロット

# %%
"""与えられたデータに対するP値函数と信頼区間のプロット"""
function plot_pvalue_funcs(n, k, α = 0.05)
    m = k/n
    s = √(m*(1 - m)/n)
    a, b = max(0, m - 4s), min(1, m + 4s)
    θ = range(a, b; length=1000)
    p_exact = pvalue_exact.(n, k, θ)
    p_dos   = pvalue_dos.(n, k, θ)
    
    ci_exact = confidence_interval(pvalue_exact, n, k, α) |> collect
    ci_dos   = confidence_interval(pvalue_dos,   n, k, α) |> collect
    
    plot(;)
    plot!(θ, p_exact; label="pvalue_exact")
    plot!(θ, p_dos; label="pvalue_dos", ls=:dash)
    plot!(ci_exact, fill(α + 0.005, 2); label="CI for pvalue_exact", c=:blue, lw=2.5)
    plot!(ci_dos,   fill(α - 0.005, 2); label="CI for pvalue_dos", c=:red, lw=2.5, ls=:dash)
    plot!(; ytick=0:0.05:1)
    plot!(; xlabel="θ", ylabel="P-value")
    title!("n = $n, k = $k, α = $α"; titlefontsize=12)
end

# %%
plot_pvalue_funcs(10, 4)

# %%
plot_pvalue_funcs(100, 40)

# %%
plot_pvalue_funcs(1000, 400)

# %% [markdown]
# ### 第一種の過誤が生じる確率

# %%
"""有効サポート"""
function effective_support(bin::Binomial; m = 4)
    μ, σ = mean(bin), std(bin)
    a = max(minimum(bin), round(Int, μ - m*σ))
    b = min(maximum(bin), round(Int, μ + m*σ))
    a:b
end

"""二項分布における期待値"""
@memoize function expectation_value(f, n, θ)
    bin = Binomial(n, θ)
    sum(k -> f(k)*pdf(bin, k), effective_support(bin))
end

"""二項分布における確率"""
probability(f, n, θ) = expectation_value(f, n, θ)

"""第一種の過誤が生じる確率をαを動かしてプロット"""
function plot_prob_typeIerror(n, θ)
    p_exact(α) = probability(k -> pvalue_exact(n, k, θ) ⪅ α, n, θ)
    p_dos(α)   = probability(k -> pvalue_dos(n, k, θ) ⪅ α,   n, θ)
    
    α = 0:0.001:1
    tick = 0:0.1:1
    P = plot(; legend=:topleft)
    plot!(α, p_exact; label="pvalue_exact")
    plot!(α, p_dos;   label="pvalue_dos", ls=:dash)
    plot!([0, maximum(α)], [0, maximum(α)]; label="", c=:black, ls=:dot)
    plot!(; xtick=tick, ytick=tick)
    plot!(; xlabel="α", ylabel="probabilty of type I error")

    α = 0:0.0001:0.1
    tick = 0:0.01:0.1
    Q = plot(; legend=:topleft)
    plot!(α, p_exact; label="pvalue_exact")
    plot!(α, p_dos;   label="pvalue_dos", ls=:dash)
    plot!([0, maximum(α)], [0, maximum(α)]; label="", c=:black, ls=:dot)
    plot!(; xtick=tick, ytick=tick)
    plot!(; xlabel="α", ylabel="probabilty of type I error")
    
    plot(P, Q; size=(800, 400), leftmargin=3Plots.mm)
end

# %% [markdown]
# αごとに第一種の過誤が実際に起こる確率のグラフが45度線に近いP値函数の方がより正確であると考えられる.
#
# 以下を見れば分かるように, `pvalue_exact` の方が `pvalue_dos` よりも正確である.

# %%
@time plot_prob_typeIerror(10, 0.4)

# %%
@time plot_prob_typeIerror(100, 0.4)

# %%
@time plot_prob_typeIerror(1000, 0.4)

# %% [markdown]
# ### 信頼区間の被覆確率

# %%
"""区間にθが含まれるかどうかを判定する函数"""
contains(interval, θ) = first(interval) ≤ θ ≤ last(interval)

"""信頼区間に真の値が含まれる確率をプロット"""
function plot_coverage_prob(n, α = 0.05, ytick=0:α/5:1)
    c_exact(θ) = probability(k -> contains(confidence_interval(pvalue_exact, n, k, α), θ), n, θ)
    c_dos(θ)   = probability(k -> contains(confidence_interval(pvalue_dos,   n, k, α), θ), n, θ)
    
    θ = 0:0.001:1
    y_exact = similar(θ)
    y_dos   = similar(θ)
    Threads.@threads for i in eachindex(θ)
        y_exact[i] = c_exact(θ[i])
        y_dos[i]   = c_dos(θ[i])
    end
    
    plot(; legend=:bottomright, ylim=(1 - 1.25α, 1 + 0.05α))
    plot!(θ, y_exact; label="pvalue_exact")
    plot!(θ, y_dos;   label="pvalue_dos", ls=:dash)
    hline!([1 - α]; label="", c=:black, ls=:dot)
    plot!(; xtick=0:0.1:1, ytick)
    plot!(; xlabel="θ", ylabel="coverage probability")
    title!("coverage probability for n = $n, α = $α"; titlefontsize=10)
end

# %% [markdown]
# 信頼区間に含まれる確率を被覆確率(coverage probability)と呼ぶ.
#
# 被覆確率が $1-\alpha$ に近い信頼区間の方がより正確であると考えられる.
#
# 以下を見れば分かるように, pvalue_exact の信頼区間の方が pvalue_dos の信頼区間(Clopper-Pearsonの信頼区間)よりも正確である.

# %%
@time plot_coverage_prob(10)

# %%
@time plot_coverage_prob(100)

# %%
@time plot_coverage_prob(1000)

# %%
