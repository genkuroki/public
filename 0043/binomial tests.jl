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
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# # 二項検定達
#
# * 黒木玄
# * 2023-06-09
# $
# \newcommand\op{\operatorname}
# \newcommand\cdf{\op{cdf}}
# \newcommand\ccdf{\op{ccdf}}
# \newcommand\confint{\op{confint}}
# \newcommand\pvalue{\op{pvalue}}
# \newcommand\Binomial{\op{Binomial}}
# \newcommand\Beta{\op{Beta}}
# \newcommand\Normal{\op{Normal}}
# \newcommand\phat{\hat{p}}
# \newcommand\quantile{\op{quantile}}
# $
#
# ## P値の4種類の定義の仕方
#
# 二項分布モデル
#
# $$
# P(k|n,p) = \binom{n}{k}p^k(1-p)^{n-k} \quad (k=0,1,2,\ldots,n)
# $$
#
# に関するパラメータ $0\le p\le 1$ に関する両側検定のP値の4種類の定義の仕方を紹介しよう(他にも定義の仕方ある). 
#
# 分布 $D$ に従う確率変数を $X$ と書く. 分布 $D$ 関する $\cdf$, $\ccdf$ を次のように定める:
#
# $$
# \cdf(D, x) = P(X\le x), \quad
# \ccdf(D, x) = 1 - \cdf(D, x) = P(X > x).
# $$
#
# ### Clopper-Pearsonの信頼区間を与えるP値
#
# 二項分布内でデータの数値 $k$ 以下になる確率とデータの数値 $k$ 以下になる確率の小さい方の2倍をP値の定義とする:
#
# $$
# \pvalue_{\text{Clopper-Pearson}}(k|n,p) = 
# \min\left(
# 1,
# 2\cdf(\Binomial(n,p), k), 
# 2\ccdf(\Binomial(n,p), k-1) 
# \right).
# $$
#
# ただし, P値の値を1以下に制限するために少しだけ改変していることに注意せよ.
#
# この定義の仕方は片側検定のP値の2倍を両側検定のP値として採用する方針になっている.
#
# ### Sterneの信頼区間を与えるP値
#
# $$
# \pvalue_{\text{Sterne}}(k|n,p) = \sum_{P(j|n,p)\le P(k|n,p)} P(j|n,p).
# $$
#
# 二項分布内で生じる確率がデータの数値 $k$ が生じる確率以下の $j$ の確率の和をP値の定義としている.
#
# この方法では自動的に両側の確率の和が計算される.
#
# ## Wilsonの信頼区間を与えるP値
#
# 中心極限定理によって, 二項分布 $\Binomial(n,p)$ はそれと同じ平均 $np$ と標準偏差 $\sqrt{np(1-p)}$ を持つ正規分布 $\Normal(np, \sqrt{np(1-p)})$ で近似される.  その近似を使って両側検定のP値を次のように定義できる:
#
# $$
# \pvalue_{\text{Wilson}}(k|n,p) = \min\left(
# 2\cdf(\Normal(np, \sqrt{np(1-p)}), k), 
# 2\ccdf(\Normal(np, \sqrt{np(1-p}), k) 
# \right).
# $$
#
# ## Waldの信頼区間を与えるP値
#
# $\phat = k/n$ とおき, WilsonのP値の定義で使った二項分布の標準偏差 $\sqrt{np(1-p)}$ をデータの数値から計算される $\sqrt{n\phat(1-\phat)}$ で置き換えて次のように両側検定のP値を定義する:
#
# $$
# \pvalue_{\text{Wald}}(k|n,p) = \min\left(
# 2\cdf(\Normal(np, \sqrt{n\phat(1-\phat)}), k), 
# 2\ccdf(\Normal(np, \sqrt{n\phat(1-\phat}), k) 
# \right).
# $$
#
# 正確な標準偏差をデータの数値から得られる標準偏差の推定値で置き換えてしまっているので, この方法で定義されたP値の性質はあまりよくない.  しかし, このP値に対応する信頼区間の計算が非常に易しくなるという利点がある.
#
# ## 信頼区間の定義
#
# P値函数を $\pvalue(k|n,p)$ と書くとき, 信頼水準 $1-\alpha$ の信頼区間は次のように定義される:
#
# $$
# \confint(k|n,\alpha) = \{\, p \in [0,1] \mid \pvalue(k|n,p) \ge \alpha\,\}.
# $$
#
# すなわち, P値が $\alpha$ 以上になるパラメータ $p$ の値全体の集合として信頼区間が一般的に定義される.
#
# ### 練習問題
#
# __問題:__ Clopper-Pearson, Wilson, Waldの場合に信頼区間の具体的な表示を求めよ. (さらに, よく知られているものに一致することを確認せよ.)
#
# __ヒント1:__ 二項分布の $\cdf$, $\ccdf$ はベータ分布の $\ccdf$, $\cdf$ で表示できる:
#
# $$
# \begin{aligned}
# &
# \ccdf(\Binomial(n, p), k-1) = \cdf(\Beta(k, n-k+1), p),
# \\ &
# \cdf(\Binomial(n, p), k) = \ccdf(\Beta(k+1, n-k), p).
# \end{aligned}
# $$
#
# 具体的な数式で書くと,
#
# $$
# \begin{aligned}
# &
# \sum_{j=k}^n \binom{n}{j} p^j(1-p)^{n-j} = \frac{1}{B(k,n-k+1)}\int_0^p t^{k-1}(1-t)^{n-k}\,dt,
# \\ &
# \sum_{j=0}^k \binom{n}{j} p^j(1-p)^{n-j} = \frac{1}{B(k+1,n-k)}\int_p^1 t^{k}(1-t)^{n-k-1}\,dt.
# \end{aligned}
# $$
#
# まず, これを証明せよ. 連続な確率分布 $D$ の累積分布函数 $p = \cdf(D, x)$ の逆函数を分位点函数と呼び $x = \quantile(D, p)$ と書く.  Clopper-Pearsonの信頼区間はベータ分布の分位点函数で表される.
#
# __ヒント1':__ ベータ分布 $\Beta(k, n-k+1)$ は一様分布 $\op{Uniform}(0,1)$ のサイズ $n$ の標本の中で $k$ 番目に小さな値(順序統計量)の分布になっている. そのことを認めれば上の二項分布とベータ分布の関係は自明になる.  一様分布のサイズ $n$ の標本の中で $k$ 番目に小さな値が $p$ 以下になることと, 一様分布のサイズ $n$ の標本中に $p$ 以下のものが $k$ 個以上存在することは同値である.  ゆえに前者の確率 $\cdf(\Beta(k, n-k+1), p)$ と後者の確率 $\ccdf(\Binomial(n,p), k-1)$ (二項分布 $\Binomial(n,p)$ に従う確率変数の値が $k$ 以上になる確率)は等しい.
#
# __ヒント2:__ Wilsonの信頼区間は二次方程式の解を使って表される. 正規分布近似を使っているので, 正規分布の分位点函数 $\quantile(\Normal(\mu, \sigma), p)$ を自由に使ってよい. 実際には $z_{\alpha/2} = \quantile(\Normal(0, 1), 1-\alpha)$ のみを使えば十分である.
#
# __ヒント3:__ 最も易しいのはWaldの信頼区間の表示である.  $z_{\alpha/2} = \quantile(\Normal(0, 1), 1-\alpha)$ を自由に使ってよい.

# %%
using Distributions
using StatsPlots
default(fmt=:png)

# %%
n, k = 10, 6
plot(p -> cdf(Beta(k, n-k+1), p), 0, 1; label="")
plot!(p -> ccdf(Binomial(n, p), k-1), 0, 1; label="", ls=:dash)

# %%
n, k = 10, 6
plot(p -> ccdf(Beta(k+1, n-k), p), 0, 1; label="")
plot!(p -> cdf(Binomial(n, p), k), 0, 1; label="", ls=:dash)

# %%
function pvalue_clopper_pearson(k, n, p)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

pvalue_clopper_pearson(6, 10, 1/4)

# %%
2cdf(Beta(6, 10-6+1), 1/4)

# %%
x ⪅ y = x < y || x ≈ y

function pvalue_sterne(k, n, p)
    bin = Binomial(n, p)
    sum(pdf(bin, j) for j in support(bin) if pdf(bin, j) ⪅ pdf(bin, k))
end

pvalue_sterne(6, 10, 1/4)

# %%
function pvalue_wilson(k, n, p)
    bin = Binomial(n, p)
    μ, σ = mean(bin), std(bin)
    normal = Normal(μ, σ)
    min(1, 2cdf(normal, k), 2ccdf(normal, k))
end

pvalue_wilson(6, 10, 1/4)

# %%
function pvalue_wald(k, n, p)
    bin = Binomial(n, p)
    p̂ = k/n
    μ, σ = mean(bin), √(n * p̂ * (1 - p̂))
    normal = Normal(μ, σ)
    min(1, 2cdf(normal, k), 2ccdf(normal, k))
end

pvalue_wald(6, 10, 1/4)

# %%
using Roots

# 手抜き
function confint(pvaluefunc, k, n; α=0.05)
    find_zeros(p -> pvaluefunc(k, n, p) - α, 0, 1)
end

for t in (:clopper_pearson, :sterne, :wilson, :wald)
    f = Symbol(:confint_, t)
    g = Symbol(:pvalue_, t)
    @eval $f(k, n; α=0.05) = confint($g, k, n; α)
end

for t in (:clopper_pearson, :sterne, :wilson, :wald)
    f = Symbol(:confint_, t)
    @eval @show $f(6, 10; α=0.05)
end

# %%
for (k, t) in enumerate((:clopper_pearson, :sterne, :wilson, :wald))
    P = Symbol(:P, k)
    g = Symbol(:pvalue_, t)
    @eval $P = plot(p -> $g(6, 10, p), 0, 1; label="")
    plot!(xtick=0:0.1:1, ytick=0:0.1:1)
    title!("$t")
end
plot(P1, P2, P3, P4; size=(800, 500))

# %%
