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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# https://x.com/reel_47/status/1810901570571386994
#
# <img src="IMG_4974.png" width=70%>

# %%
using Distributions
using HypothesisTests

Red = [18, 22, 81, 45, 14, 31, 37, 183, 401]
Gold = [11, 73, 14, 18, 71, 25, 15, 2, 77, 24, 114, 26, 17, 11, 159, 25, 47, 78, 95]

@show mean(Red), std(Red)
@show mean(Gold), std(Gold)
UnequalVarianceTTest(Red, Gold)

# %%
ExactMannWhitneyUTest(Red, Gold)

# %%
ApproximateMannWhitneyUTest(Red, Gold)

# %%
function brunner_munzel_test(X, Y; p=1/2)
    m, n = length(X), length(Y)
    phat = mean((x < y) + (x == y)/2 for x in X, y in Y)
    Hbarx = n*(1 - phat)
    Hbary = m*phat
    sx2 = 1/n^2 * 1/(m-1) * sum(x -> (sum((y < x) + (y == x)/2 for y in Y) - Hbarx)^2, X)
    sy2 = 1/m^2 * 1/(n-1) * sum(y -> (sum((x < y) + (x == y)/2 for x in X) - Hbary)^2, Y)
    sehat = √(sx2/m + sy2/n)
    tvalue = (phat - p)/sehat
    df = (sx2/m + sy2/n)^2 / ((sx2/m)^2/(m-1) + (sy2/n)^2/(n-1))
    pvalue = sehat > 0 ? 2ccdf(TDist(df), abs(tvalue)) : phat ≈ p ? 1.0 : 0.0
    (; phat, sehat, tvalue, df, pvalue, p)
end

brunner_munzel_test(Red, Gold) |> pairs |> Dict

# %%
