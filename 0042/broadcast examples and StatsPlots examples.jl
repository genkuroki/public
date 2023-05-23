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
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
using Plots

n = 100
x = range(-1, 1, n)
u = 0.5randn(n)
y = sinpi.(x) + u

d = 3
A = x .^ (0:d)'
β̂ = A \ y
f̂(x) = evalpoly(x, β̂)

scatter(x, y; label="data")
plot!(f̂; label="degree-$d polynomial")

# %%
(-5:5) .^ (0:3)'

# %%
using Distributions
using StatsPlots

n = 10^6
X = sort!(randn(n))
Y = X + reverse(X)
@show normal = fit(Normal, Y)
stephist(Y; norm=true, label="data")
plot!(normal; label="normal approx.")

# %%
using Distributions
using StatsPlots

n = 10^6
X = sort!(randn(n))
@views Y = X[1:2:end] + X[2:2:end]
@show normal = fit(Normal, Y)
stephist(Y; norm=true, label="data")
plot!(normal; label="normal approx.")

# %%
using Distributions
using StatsPlots

n = 10^6
X = sort!(randn(n))
@views Y = X[1:end÷2] + X[end÷2+1:end]
@show normal = fit(Normal, Y)
stephist(Y; norm=true, label="data")
plot!(normal; label="normal approx.")

# %%
