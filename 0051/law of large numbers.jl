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
#using StatsPlots
using Plots
default(fmt=:png)

n = 10^4
L = 100
X = rand(TDist(2), n, L)
Y = cumsum(X; dims=1) ./ (1:n)
ns = [1:10; n÷100:n÷100:n]
plot(ns, Y[ns,:]; label="", lw=0.3, alpha=0.6, ylim=(-1, 1))
plot!(n -> 6/√n; label="±6/√n", c=:red)
plot!(n -> -6/√n; label="", c=:red)

# %%
