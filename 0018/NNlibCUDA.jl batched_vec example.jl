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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
using CUDA
using StaticArrays

# %%
θ = range(0, 2π; length=361)[1:end-1]
θ_cu = cu(collect(θ))

# %%
r_smf32(t) = @SMatrix(Float32[
    cos(t) -sin(t) 0
    sin(t)  cos(t) 0
    0      0       1
])

f(t) = r_smf32(t) * @SVector Float32[1, 0, 0]

# %%
f.(θ_cu)

# %%
using CUDA, NNlib, NNlibCUDA

# %%
?batched_vec

# %%
R(t) = [
    cos(t) -sin(t) 0
    sin(t)  cos(t) 0
    0      0       1
]

θ = range(0, 2π; length=361)[1:end-1]
A = mapslices(t -> R(t[1, 1]), reshape(θ, 1, 1, :); dims=(1, 2))
A_cu = cu(A)

# %%
batched_vec(A_cu, cu([1.0, 0.0, 0.0]))

# %%
using CUDA, NNlib, NNlibCUDA

R(t) = [
    cos(t) -sin(t) 0
    sin(t)  cos(t) 0
    0      0       1
]

θ = range(0, 2π; length=361)[1:end-1]
A = mapslices(t -> R(t[1, 1]), reshape(θ, 1, 1, :); dims=(1, 2))
batched_vec(cu(A), cu([1.0, 0.0, 0.0]))

# %% [markdown]
# See also https://discourse.julialang.org/t/how-to-broadcast-or-batch-multiply-a-batch-of-matrices-with-another-matrix-on-the-gpu/67259

# %%
