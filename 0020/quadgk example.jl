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
using QuadGK

f(x) = (x + 2) * x * (x - 2)
g(x) = -f(5 - x)

@show quadgk(f, 0, 5)
@show quadgk(f, 5, 0)
@show quadgk(g, 0, 5);

# %%
