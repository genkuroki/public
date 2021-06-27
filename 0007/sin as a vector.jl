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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
module O

Base.@kwdef struct Foo
    d::Dict{Float64, Float64} = Dict{Float64, Float64}()
end

Base.getindex(foo::Foo, x::Float64) = get(foo.d, x, sin(x))
Base.getindex(foo::Foo, x::Int) = getindex(foo, float(x))
Base.getindex(foo::Foo, x::AbstractArray) = getindex.(Ref(foo), x)

Base.setindex!(foo::Foo, v, x::Float64) = setindex!(foo.d, v, x)
Base.setindex!(foo::Foo, v, x::Int) = setindex!(foo, v, float(x))
Base.setindex!(foo::Foo, v, x::AbstractArray) = setindex!.(Ref(foo), v, x)

end

# %%
foo = O.Foo()

# %%
foo[π/6]

# %%
foo[π/6] = 30
foo[π/6]

# %%
foo[1]

# %%
foo[1] = -99
foo[1]

# %%
foo[range(0, π/2; length=7)]

# %%
foo[range(0, π/2; length=7)] = 0:6
foo[range(0, π/2; length=7)]

# %%
foo[π/6]

# %%
