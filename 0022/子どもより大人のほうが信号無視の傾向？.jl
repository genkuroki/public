# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# * https://www3.nhk.or.jp/kansai-news/20210921/2000051589.html
# * [public/0014/2×2の分割表の独立性検定の信頼区間.ipynb](https://github.com/genkuroki/public/blob/main/0014/2%C3%972%E3%81%AE%E5%88%86%E5%89%B2%E8%A1%A8%E3%81%AE%E7%8B%AC%E7%AB%8B%E6%80%A7%E6%A4%9C%E5%AE%9A%E3%81%AE%E4%BF%A1%E9%A0%BC%E5%8C%BA%E9%96%93.ipynb)
#
# ![2021-09-23.png](2021-09-23.png)

# %%
using Distributions

A = [
     14-0   0
    833-78 78
    104-12 12
]

p = sum(A[:,2])/sum(A)

# %%
null = product_distribution([Binomial(sum(A[i,:]), p) for i in axes(A, 1)])

# %%
L = 10^6
X = rand(null, L)

# %%
a = mapslices(v -> v[1] ≤ A[1,2] && v[2]+v[3] ≥ A[2,2]+A[3,2], X; dims=1)

# %%
mean(a)

# %% [markdown]
# 仮に子供も大人もお年寄りも信号無視する割合が同じ9.5%程度であったときに、子供14人、大人833人、お年寄り104人の調査で、子供が誰も信号無視せず、大人とお年寄りで信号無視する人達の人数の合計が78人+12人以上になる確率は、11％程度になる。

# %%
A = [
     14-1   1
    833-78 78
    104-12 12
]
q = sum(A[:,2])/sum(A)

null2 = product_distribution([Binomial(sum(A[i,:]), p) for i in axes(A, 1)])
L = 10^6
X2 = rand(null2, L)
a2 = mapslices(v -> v[1] ≤ A[1,2] && v[2]+v[3] ≥ A[2,2]+A[3,2], X; dims=1)
mean(a2)

# %% [markdown]
# 上と同じ計算を、仮に子供14人中信号無視をしてしまった子が１人出てしまった場合に行ってみると、結果は28%程度になる。

# %%
