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
mutable struct TypeA
    n::Int
end

mutable struct TypeB
    x::Float64
end

"""
`fun!(a::TypeA, b::TypeB, c::AbstractVector{Float64})` only mutates `b`.
"""
function fun!(a::TypeA, b::TypeB, c::AbstractVector{Float64})
    b.x = sum(c) / a.n
end

a, b, c = TypeA(10), TypeB(0.0), collect(1.0:10.0)
fun!(a, b, c)
@show a b c;

# %%
@doc fun!

# %%
?fun!

# %%
