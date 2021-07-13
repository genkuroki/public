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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using Plots
pyplot()

function halftorus(n=25, a=5, b=10; lims=(-11, 11), size=(500, 400))
    u = range(-π/2, π/2, length=n+1)
    v = range(0, 2π, length=2n+1)
    x = @. (b + a * cos(u')) * cos(v)
    y = @. (b + a * cos(u')) * sin(v)
    z = @.      a * sin(u')  * one(v)
    surface(x, y, z; lims, colorbar=false, size)
end

halftorus()

# %%
