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
#     display_name: Julia 1.7.0-beta2
#     language: julia
#     name: julia-1.7
# ---

# %%
@show VERSION

using Plots

f(x, y) = (1 - x)^2 + 100(y - x^2)^2

X = Y = range(-5, 5; length=1001)
@show x, y = argmin(Base.splat(f), Iterators.product(X, Y))
@show f(x, y) 

surface(X, Y, log∘f; color=:rainbow, clim = (-7.5, 10), camera=(30, 60))
scatter!([x], [y], [-7.5]; label="minimum", color=:cyan)
plot!(; xtick = -5:5, ytick = -5:5, zlim=(-7.5, 10))

# %%
