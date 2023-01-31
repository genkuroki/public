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
#     display_name: Julia 1.9.0-beta3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png,
    titlefontsize=8, tickfontsize=6, guidefontsize=7, legendfontsize=7,
    plot_titlefontsize=8,
    size=(400, 250))

# 離散分布のグラフを書くための設定
mypdf(dist, x) = pdf(dist, x)
mypdf(dist::DiscreteUnivariateDistribution, x) = pdf(dist, round(Int, x))

# 分布名から {～} を削除する
distname(dist) = replace(string(dist), r"{[^}]*}"=>"")

# 尖度
mykurtosis(dist) = kurtosis(dist)

# %%
Normal()

# %%
distname(Normal())

# %%
# 分布の設定
dist = Gamma(2, 3)

# %%
# 分布のプロット
xlim = (-0.5, 30)
plot(x -> pdf(dist, x), xlim...; label="", title=distname(dist))

# %%
# 分布のサンプルを大量に生成
n = 100
L = 10^5
samples = [rand(dist, n) for _ in 1:L]
first(samples, 3)

# %%
# サンプル達の標本平均と不偏分散を計算してヒストグラムをプロット
means_of_samples = mean.(samples)
variances_of_samples = var.(samples)

println("(dist, n) = (", distname(dist), ", ", n, ")")

bin = 100

P1 = stephist(means_of_samples; norm=true, bin, label="")
title!("means of samples")
plot!(Normal(mean(dist), std(dist)/√n); label="CLT", ls=:dash)

P2 = stephist((n-1)*variances_of_samples/var(dist); norm=true, bin, label="")
title!("(n-1)×(variances of samples)/var(dist)")
ku = kurtosis(dist)
plot!(Normal(n-1, (n-1)*√(ku/n + 2/(n-1))); label="CLT", ls=:dash)
plot!(Chisq(n-1); label="Chisq(n-1)", ls=:dashdot)

plot(P1, P2; size=(800, 250))

# %% [markdown]
# $X_1,X_2,\ldots,X_n$ は分布 $D$ のi.i.d.であるとし,
#
# $$
# \bar{X} = \frac{1}{n}\sum_{i=1}^n X_i, \quad
# S^2 = \frac{1}{n-1}\sum_{i=1}^n (X_i - \bar{X})^2
# $$
#
# とおく. このとき, 中心極限定理より, $n$ を十分に大きくすると,
#
# * $\bar{X}$ の分布は近似的に $\operatorname{Normal}(\mu, \sigma/\sqrt{n})$ になる.
# * $(n-1)S^2/\sigma^2$ の分布は近似的に $\operatorname{Normal}\left(n-1, (n-1)\sqrt{\bar{\kappa}_4/n + 2/(n-1)}\,\right)$ になる.
#
# ただし, ここで $\mu$ と $\sigma^2$ は分布 $D$ の期待値と分散であり, $\bar{\kappa}_4$ は分布 $D$ の(過剰)尖度である.
#
# さらにもしも $D$ が正規分布ならば,
#
# * $\bar{X}$ の分布はぴったり $\operatorname{Normal}(\mu, \sigma/\sqrt{n})$ になる.
# * $(n-1)S^2/\sigma^2$ の分布はぴったり $\operatorname{Chisq}(n-1)$ になる.
#
# 正規分布に関する後者の結果は分散の値に関する検定や信頼区間で使われる.
#
# 分布 $\operatorname{Chisq}(n-1)$ の期待値は $n-1$ であり, 分布 $\operatorname{Normal}\left(n-1, (n-1)\sqrt{\bar{\kappa}_4/n + 2/(n-1)}\,\right)$ の期待値 $n-1$ に一致する.
#
# しかし, 分布 $\operatorname{Chisq}(n-1)$ の分散は $2\sqrt{n-1}$ であり, 分布 $\operatorname{Normal}\left(n-1, (n-1)\sqrt{\bar{\kappa}_4/n + 2/(n-1)}\,\right)$ の分散 $\bar{\kappa}_4(n-1)/n + 2(n-1)$ との差は $n$ が大きなときほぼ分布 $D$ の過剰尖度 $\bar{\kappa}_4$ に等しくなる.
#
# ゆえに, $(n-1)S^2/\sigma^2$ の分布はぴったり $\operatorname{Chisq}(n-1)$ になることを使う分散に関する検定や信頼区間は, たとえどんなに $n$ を大きくしても, 分布 $D$ の過剰尖度 $\bar{\kappa}_4$ の分だけ不正確になる.
#
# それとは対照的に, $\bar{X}$ の分布はぴったり $\operatorname{Normal}(\mu, \sigma/\sqrt{n})$ になることを使う平均に関する検定や信頼区間は $n$ を十分に大きくすると中心極限定理によって誤差が小さくなってくれる.
#
# このように, 正規分布モデルを使った検定や信頼区間であっても, 平均を扱う場合と分散を扱う場合では頑健さが大違いになる.

# %%
# 以上を函数化

function plot_means_and_vars(dist, n, L=10^5; bin=100, vbin=bin, vlim=nothing)
    samples = [rand(dist, n) for _ in 1:L]
    
    means_of_samples = mean.(samples)
    variances_of_samples = var.(samples)

    P1 = stephist(means_of_samples; norm=true, bin, label="")
    title!("means of samples")
    plot!(Normal(mean(dist), std(dist)/√n); label="CLT", ls=:dash)

    P2 = plot()
    if isnothing(vlim)
        stephist!((n-1)*variances_of_samples/var(dist); norm=true, bin=vbin, label="")
    else
        stephist!((n-1)*variances_of_samples/var(dist); norm=true, bin=vbin, label="", xlim=vlim)
    end
    title!("(n-1)×(variances of samples)/var(dist)")
    ku = kurtosis(dist)
    plot!(Normal(n-1, (n-1)*√(ku/n + 2/(n-1))); label="CLT", ls=:dash)
    plot!(Chisq(n-1); label="Chisq(n-1)", ls=:dashdot)
    
    plot(P1, P2; size=(800, 250))
    plot!(plot_title="$(distname(dist)), n=$n, kurtosis=$(round(ku; digits=3))")
end

# %%
plot_means_and_vars(Normal(2, 3), 20)

# %%
plot_means_and_vars(Normal(2, 3), 100)

# %%
plot_means_and_vars(Gamma(2, 3), 20)

# %%
plot_means_and_vars(Gamma(2, 3), 100)

# %%
plot_means_and_vars(Gamma(10, 3), 100)

# %%
plot_means_and_vars(Gamma(100, 3), 100)

# %%
plot(Gamma(100, 3), 180, 500; label=distname(Gamma(100, 3)))
plot!(Normal(mean(Gamma(100, 3)), std(Gamma(100, 3))); label="normal approx.", ls=:dash)
#plot!(legend=:topleft)

# %%
plot_means_and_vars(Exponential(), 100)

# %%
plot_means_and_vars(Laplace(), 100)

# %%
plot_means_and_vars(TDist(4.1), 100; vbin=1000, vlim=(-200, 400))

# %%
plot_means_and_vars(TDist(4.1), 1000; vbin=1000, vlim=(0, 2000))

# %%
plot_means_and_vars(TDist(5), 100; vlim=(0, 200), vbin=1000)

# %%
plot_means_and_vars(TDist(10), 100; vlim=(0, 200), vbin=100)

# %%
plot_means_and_vars(Uniform(), 100)

# %%
plot_means_and_vars(Bernoulli(0.3), 100; bin=20)

# %%
plot_means_and_vars(Beta(0.3, 0.7), 100)

# %%
plot(Beta(0.3, 0.7); ylim=(-0.1, 5.1), label=distname(Beta(0.1, 0.2)))
plot!(legend=:top)

# %%
