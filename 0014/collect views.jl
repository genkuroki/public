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

# %% [markdown]
# https://docs.julialang.org/en/v1/base/arrays/#Base.PermutedDimsArrays.PermutedDimsArray

# %%
using ArraysOfArrays
f(inp) = nestedview(PermutedDimsArray(inp, (2, 3, 1)), 2)
w, h = 5, 5
inp = 0.1(1:3) .+ 10(1:5)' .+ reshape(1:5, 1, 1, :)
out = f(inp)

# %%
out[2]

# %%
using BenchmarkTools
w, h = 1000, 1000
inp = 0.1(1:3) .+ 10(1:5)' .+ reshape(1:5, 1, 1, :)
@show [@view(inp[i, :, :]) for i in 1:3] == collect(eachslice(inp; dims=1)) == f(inp)
@btime [@view($inp[i, :, :]) for i in 1:3]
@btime collect(eachslice($inp; dims=1))
@btime f($inp)
@btime PermutedDimsArray($inp, (2, 3, 1));

# %%
