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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %%
using OffsetArrays

# %%
v = OffsetVector(rand(0:9, 9), -4:4)
axes(v, 1)

# %%
enumerate(v) |> collect

# %%
pairs(v)

# %%
keys(v)

# %%
A = OffsetArray(rand(0:9, 5, 5), 0:4, -2:2)
axes(A)

# %%
enumerate(A) |> collect

# %%
pairs(A) |> collect

# %%
keys(A)

# %%
