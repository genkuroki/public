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
#     display_name: Julia 1.11.1
#     language: julia
#     name: julia-1.11
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

apply(f, x...) = f(x...)

# %%
dist(μ) = μ + TDist(3)

#n = 20
#X = rand(dist(0), n)
@show mean(X)
f(μ) = exp(loglikelihood(dist(μ), X) - loglikelihood(dist(0), X))
plot(f, -3, 3)
vline!([mean(X)])

# %%
dist(μ) = μ + TDist(3)

n = 10
X = rand(dist(0), n)
g(x) = exp(loglikelihood(dist(1), x .+ X) - loglikelihood(dist(0), x .+ X))
plot(g, -5, 25)

# %%
X

# %%
@show X;

# %%
histogram(X)

# %%
@show round.(X; digits=1);

# %%
Y = [1.88, -0.361, 0.594, 2.19, -0.364, -2.12, 0.638, -0.503, -6.48, -0.158]

@show mean(Y)
f(μ) = exp(loglikelihood(dist(μ), Y) - loglikelihood(dist(0), Y))
plot(f, -3, 3)
vline!([mean(Y)])

# %%
dotplot(Y; size=(200, 400), xtick=1:1, xlim=(0.5, 1.5))

# %%
