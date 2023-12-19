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
#     display_name: Julia 1.9.4
#     language: julia
#     name: julia-1.9
# ---

# %%
# Julia言語

using Distributions
using StatsPlots
default(fmt=:png)

dist = Gamma(2, 3)
n = 10000
X = rand(dist, n)
μ = mean(X)
σ = std(X)
Z = @. (X - μ)/σ

binx = 0:0.5:40
binz = @.(binx - μ)/σ
P = histogram(X; label="original data X", bin=binx, c=1)
Q = histogram(Z; label="Z-score normalization of X", bin=binz, c=2)
plot(P, Q; layout=(2, 1))

# %%
