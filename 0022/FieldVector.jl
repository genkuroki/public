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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using StaticArrays
struct Point{T} <: FieldVector{2, T} x::T; y::T end

# %%
p = Point(1, 2)

# %%
p[2], p.y

# %%
q = Point(3, 4)

# %%
q - p

# %%
3p

# %%
using LinearAlgebra
dot(p, q), cross(p, q)

# %%
struct PointPlain{T} x::T; y::T end
methodswith(PointPlain; supertypes=true)

# %%
struct PointFieldVector{T} <: FieldVector{2, T} x::T; y::T end
methodswith(PointFieldVector; supertypes=true)

# %%
