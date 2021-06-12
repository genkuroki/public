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

# %% tags=[]
surface(x, y, f; xlim=extrema(x), ylim=extrema(y), camera=(30, 45), size=(500, 500), colorbar=false)
scatter!([pt[1]], [pt[2]], [f(pt...)]; label="", color=:cyan, ms=2)

# %%
f(p) = sin(p[1]) * sin(p[2]) + (p[1] - 2p[2])^2/50
x = y = range(-π, π; length=401)
findmin(f, Iterators.product(x, y))

# %%
f(p) = sin(p[1]) * sin(p[2]) + (p[1] - 2p[2])^2/50
x = y = range(-π, π, 1000)
findmin(f, Iterators.product(x, y))

# %%
using Plots
plotly()
x = y = range(-3, 3, length=100)
f(x, y) = exp(-x^2-y^2)
z = f.(x', y)
surface(x, y, z; colorbar=false, size=(600, 600), color=:jet)

# %%
using Plots
plotly()
x = y = range(-3, 3, length=100)
z = exp.(.-x'.^ 2 .- y.^ 2)
surface(x, y, z; colorbar=false, size=(600, 600), color=:jet)

# %%
using Plots
plotly()
x = y = range(-3, 3, length=100)
z = @. exp(-x'^2 .- y^2)
surface(x, y, z; colorbar=false, size=(600, 600), color=:jet)

# %%
