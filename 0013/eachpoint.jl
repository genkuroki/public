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
# https://discourse.julialang.org/t/generically-iterating-over-arrays-of-points/65273/2

# %%
using StaticArrays

eachpoint(itr) = itr
eachpoint(ptmatrix::AbstractMatrix) = eachcol(ptmatrix)

function f(pts)
    for pt in eachpoint(pts)
        println(pt)
    end
end

f(SVector(k+1, k+2, k+3) for k in 0:3:12)

# %%
f(reshape(1:15, 3, 5))

# %%
