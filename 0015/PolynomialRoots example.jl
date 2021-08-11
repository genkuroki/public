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
using PolynomialRoots
p = [1; -1; fill(0, 98); 1]
a = roots(p)
@show [(k, round(real(sum(a .^ k)); digits=6)) for k in 1:100];

# %%
