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
abstract type AbstractSep{T} end

struct Sep{T<:Union{Nothing, Tuple{<:AbstractSep, <:AbstractSep}}} <: AbstractSep{T}
    a::T
end
Sep(x::AbstractSep, y::AbstractSep) = Sep((x, y))

const Nil = Sep(nothing)
Base.show(io::IO, ::Sep{Nothing}) = print(io, "Nil")
Base.show(io::IO, x::Sep{<:Tuple}) = print(io, x.a)

# %%
Nil

# %%
a = Sep(Nil, Nil)

# %%
b = Sep(a, a)

# %%
c = Sep(a, b)

# %%
typeof(c)

# %%
