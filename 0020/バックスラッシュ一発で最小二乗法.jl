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

# %%
using Plots

f(x) = 1 - 2x
x = range(0, 2; length=20)
y = @. f(x) + randn()
X = x .^ (0:1)'

β̂ = X \ y
g(x) = β̂[1] + β̂[2]*x

xs = range(extrema(x)...; length=300)
scatter(x, y; label="data", color=1)
plot!(xs, f; label="true line", color=1, ls=:dash)
plot!(xs, g; label="regression line", color=2, lw=2)

# %%

# %%
