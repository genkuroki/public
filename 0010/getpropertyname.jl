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
struct Param{T} _inner::Dict{Symbol, T} end
Param() = Param(Dict{Symbol, Any}())
Base.getproperty(param::Param, name::Symbol) = getfield(param, :_inner)[name]
Base.setproperty!(param::Param, name::Symbol, val) = setindex!(getfield(param, :_inner), val, name)
Base.propertynames(param::Param) = (keys(getfield(param, :_inner))...,)

# %%
p = Param()

# %%
propertynames(p)

# %%
p.a = 123
p.b = 4.56
p.c = "hoge"
p

# %%
propertynames(p)

# %%
p.

# %%
