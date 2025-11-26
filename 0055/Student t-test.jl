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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
ENV["LINES"] = 200
ENV["COLUMNS"] = 200
using Distributions
using DataFrames
using HypothesisTests
using RCall

rd(x) = round(x; sigdigits=2)

dfs = 4:13
ts = 0:0.1:6
tbl_pvals = @. rd(2ccdf(TDist(dfs'), ts))
names = ["t-value", ("df=$k" for k in dfs)...]
data_pvals = DataFrame([ts tbl_pvals], names)
print(data_pvals)

# %% [markdown]
# 以上はt分布での両側P値の表。
#
# 以下の表は中原治『基礎から学ぶ統計学』p.169より
#
# <img width=350 src="IMG_1453.jpeg">

# %%
x_A = [118, 132, 120, 115, 113]
x_B = [129, 126, 134, 135, 131]
n_A = length(x_A)
n_B = length(x_B)
@show xbar_A = mean(x_A)
@show xbar_B = mean(x_B)
@show s_A = std(x_A)
@show s_B = std(x_B)
@show df = n_A + n_B - 2
@show s_p = sqrt(((n_A - 1) * s_A^2 + (n_B - 1) * s_B^2) / df)
@show t = (xbar_A - xbar_B) / (s_p * sqrt(1/n_A + 1/n_B))
;

# %% [markdown]
# ここまでは電卓で計算できる。
#
# 上の両側P値の表から df = 8 と |t| ≈ 3.1 の値を読み取ると P値 ≈ 0.015 だと分かる。
#
# 以下は表を使わない計算。

# %% [markdown]
# 次のセルはRのptによる計算

# %%
@rput df t
@show rcopy(R"""2*pt(abs(t), df, lower.tail=F)""");

# %% [markdown]
# 次のセルはRのt.testによる計算

# %%
@rput x_A x_B
R"""
t.test(x_A, x_B, var.equal=T)
"""

# %%
@rput x_A x_B
R"""
t.test(x_A, x_B)
"""

# %% [markdown]
# 以下はJuliaでの計算

# %%
@show 2ccdf(TDist(df), abs(t));

# %%
EqualVarianceTTest(x_A, x_B)

# %%
xxA = copy(x_A)
xxA[3] = 130
@show x_A xxA x_B
EqualVarianceTTest(xxA, x_B)

# %%
UnequalVarianceTTest(x_A, x_B)

# %%
xxA = copy(x_A)
xxA[3] = 130
@show x_A xxA x_B
UnequalVarianceTTest(xxA, x_B)

# %%
