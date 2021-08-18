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
#     display_name: Julia 1.7.0-beta3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Plots
x = range(-5, 5; length=300)
y = range(-3, 3; length=200)
f(x, y) = exp(-x^2/2-y^2)
z = f.(x', y)
surface(x, y, z; camera=(30, 60), color=:gist_earth)

# %%
surface(x, y, f; camera=(30, 60), color=:gist_earth)

# %%
