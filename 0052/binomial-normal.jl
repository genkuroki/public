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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# %% [markdown]
# # 二項分布と正規分布の関係
#
# $
# \newcommand\ds{\displaystyle}
# \newcommand\op{\operatorname}
# \newcommand\Binomial{\op{Binomial}}
# \newcommand\Normal{\op{Normal}}
# \newcommand\phat{\hat{p}}
# $

# %%
using Distributions
using StatsPlots
default(fmt=:png, 
    size=(600, 340), linewidth=2, 
    tickfontsize=8, guidefontsize=12, legendfontsize=12, titlefontsize=16)

# %% [markdown]
# ## 二項分布
#
# あたりが出る確率が$p$のルーレットを$n$回まわしたときにちょうど$k$回あたりが出る確率は
# $$
# {}_nC_k p^k (1 - p)^{n-k}
# $$
# になる. このとき、$k$は二項分布$\Binomial(n, p)$にしたがうという. 
#
# 例えば, $n=20$, $p=0.3$, $k=3$のとき
# $$
# {}_nC_k p^k (1 - p)^{n-k}
# = \frac{20\cdot 19\cdot 18}{3\cdot2\cdot1} 0.3^3 0.7^{17}
# = 20\cdot19\cdot3\cdot0.3^3 0.7^{17}
# \approx 0.072
# $$

# %%
20*19*3*0.3^3*0.7^17

# %%
n, p = 20, 0.3
bar(Binomial(n, p); la=0.7, fa=0.3, label="Binomial(n, p)")
#plot!(Normal(n*p, √(n*p*(1-p))); label="Normal(np, √(np(1-p))")
plot!(xguide="k", xtick=0:20, ytick=0:0.01:0.2)
title!("n = $n,  p = $p")

# %% [markdown]
# 二項分布$\Binomial(n, p)$にしたがう$k$の期待値は$\mu=np$になる. 
#
# 例えば, あたりが出る確率が$p=0.3$のルーレットを$n=20$回まわしたときに出るあたりの回数$k$の期待値は$\mu=20\cdot0.3=6$回になる.
#
# 二項分布$\Binomial(n, p)$にしたがう$k$の標準偏差(ばらつきの幅の大きさ)は$\sigma=\sqrt{np(1-p)}$になる. ($k$の分布は$\mu\pm2\sigma$のあいだに95%程度以上含まれる.)
#
# 例えば, あたりが出る確率が$p=0.3$のルーレットを$n=20$回まわしたときに出るあたりの回数$k$の標準偏差(ばらつきの幅の大きさ)は$\sigma=\sqrt{20\cdot0.3\cdot0.7}=\sqrt{4.2}\approx 2$になる. そのとき$k$の分布は2以上10以下に97.5%程度含まれる.

# %%
√(20*0.3*0.7)

# %%
cdf(Binomial(20, 0.3), 10)- cdf(Binomial(20, 0.3), 1)

# %%
n, p = 20, 0.3
bar(0:1, Binomial(n, p); c=1, la=0.7, fa=0.3, label="")
bar!(11:20, Binomial(n, p); c=1, la=0.7, fa=0.3, label="")
bar!(2:10, Binomial(n, p); c=:red, la=1, fa=0.3, label="2 ≤ k ≤ 10")
#plot!(Normal(n*p, √(n*p*(1-p))); c=2, label="Normal(np, √(np(1-p))")
plot!(xguide="k", xtick=0:20, ytick=0:0.01:0.2)
title!("n = $n,  p = $p")

# %% [markdown]
# 二項分布$\Binomial(n, p)$にしたがう$k$について$\phat=k/n$ (あたりが確率$p$で出るルーレットを$n$回まわしたときに出たあたりの回数の割合)の期待値は$\ds \frac{np}{n}=p$になり, 標準偏差は$\ds\frac{\sqrt{np(1-p)}}{n}=\sqrt{\frac{p(1-p)}{n}}$になる.

# %% [markdown]
# ## 正規分布との近似関係
#
# 期待値$\mu$と標準偏差$\sigma$を持つ正規分布$\Normal(\mu, \sigma)$の確率密度関数は
# $$
# f_\text{normal}(x|\mu,\sigma) =
# \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right),
# \quad \exp(X)=e^X.
# $$
#
# $np$と$n(1-p)$が十分大きいならば, 二項分布での確率は二項分布と同じ期待値と標準偏差を持つ正規分布の確率密度で近似されることが知られている:
# $$
# {}_nC_k p^k (1 - p)^{n-k} \approx
# f_\text{normal}(k|np,\sqrt{np(1-p)}) =
# \frac{1}{\sqrt{2\pi np(1-p)}} \exp\left(-\frac{(k-\mu)^2}{2np(1-p)}\right).
# $$
# 例えば$n=20$, $p=0.3$の場合については以下のグラフを参照せよ.

# %%
n, p = 20, 0.3
bar(Binomial(n, p); la=0.7, fa=0.3, label="Binomial(n, p)")
plot!(Normal(n*p, √(n*p*(1-p))); label="Normal(np, √(np(1-p))")
plot!(xguide="k", xtick=0:20, ytick=0:0.01:0.2)
title!("n = $n,  p = $p")

# %%
n, p = 20, 0.3
bin = Binomial(n, p)
normal = Normal(mean(bin), std(bin))
stack(Any[k, pdf(bin, k), pdf(normal, k), (pdf(normal, k) - pdf(bin, k))] for k in 0:20)'

# %%
