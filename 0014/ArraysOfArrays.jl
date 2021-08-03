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
using BenchmarkTools, ArraysOfArrays

inp = rand(5, 5, 3)
@show nestedview(inp, 2) == [inp[:, :, k] for k in axes(inp, 3)]
print("Comprehension:          ")
@btime out = [$inp[:, :, k] for k in axes($inp, 3)]
print("Comprehension with @view:")
@btime out = [@view($inp[:, :, k]) for k in axes($inp, 3)]
print("ArraysOfArrays.nestedview:")
@btime out = nestedview($inp, 2);

# %%
