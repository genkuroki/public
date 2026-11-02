# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# ## オリジナル
#
# * https://saiensu.co.jp/search/?isbn=978-4-7819-1678-1&y=2026#support
#   * code19_chap6_reidai_confidence_interval_code_Julia.txt

# %%
# パッケージがインストール済なら次の行は不要
# import Pkg; Pkg.add(["Distributions"])
using Distributions
n = 300
pvar = 210 / 300
# 二項検定の実行結果から信頼区間の項目を取得
I = quantile(Binomial(n, pvar), [0.025, 0.975])
println("pの95%信頼区間：", [I[1]/n, I[2]/n])

# %% [markdown]
# これだと「二項検定の実行結果から信頼区間の項目を取得」になっていない。
#
# ブートストラップ法で求めた区間推定の計算になってしまっている。

# %% [markdown]
# ## 訂正案

# %%
# パッケージがインストール済なら次の行は不要
# import Pkg; Pkg.add(["HypothesisTests"])
using HypothesisTests
n = 300
x = 210
level = 0.95
method = :clopper_pearson
# 二項検定の実行結果から信頼区間の項目を取得
result = BinomialTest(x, n)
println("pの$(round(Int, 100level))%信頼区間：", confint(result; level))

# %% [markdown]
# ## さらなる訂正案
#
# 普通はP値にも触れる。

# %%
# パッケージがインストール済なら次の行は不要
# import Pkg; Pkg.add(["HypothesisTests"])
using HypothesisTests
n = 300
x = 210
level = 0.95
p_test = 0.6446804955369393
# 二項検定の実行結果から信頼区間の項目を取得
result = BinomialTest(x, n, p_test)
println("pの$(round(Int, 100level))%信頼区間：", confint(result; level))
println("仮説p=$(p_test)の両側P値：", pvalue(result))

# %%
?BinomialTest

# %%
?confint(::BinomialTest)

# %%
?pvalue

# %%
