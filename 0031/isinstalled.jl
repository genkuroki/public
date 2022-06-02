# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Pkg
isinstalled(x::AbstractString) = x ∈ keys(Pkg.project().dependencies)

# %%
isinstalled("IJulia")

# %%
isinstalled("FooBarPackage")

# %%
ENV["LINES"] = 100
Pkg.project().dependencies

# %%

# %%
