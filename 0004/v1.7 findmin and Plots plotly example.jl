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
#     display_name: Julia 1.7.0-DEV
#     language: julia
#     name: julia-1.7
# ---

# %%
using Plots
plotly()
f(x, y) = sin(x) * sin(y) + (x - 2y)^2/50
x = range(-π, π; length=401)
y = range(-π, π; length=401)
surface(x, y, f; xlim=extrema(x), ylim=extrema(y), camera=(30, 45), size=(500, 500), colorbar=false)

# %%
fmin, pt = findmin(p -> f(p...), Iterators.product(x, y))

# %%
surface(x, y, f; xlim=extrema(x), ylim=extrema(y), camera=(30, 45), size=(500, 500), colorbar=false)
scatter!([pt[1]], [pt[2]], [f(pt...)]; label="", color=:cyan, ms=2)

# %%
