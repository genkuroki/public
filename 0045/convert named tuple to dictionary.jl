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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Statistics

X = randn(10)
X̄ = mean(X)
S² = std(X)
t = (; X, X̄, S²)

# %%
t.S², t[:S²], t[3]

# %%
p = pairs(t)

# %%
d = Dict(pairs(t))

# %%
