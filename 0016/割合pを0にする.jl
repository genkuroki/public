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
A = collect(reshape(1:24, 4, 6))

# %%
f(x, p) = rand() ≤ p ? zero(x) : x
B = f.(A, 0.25)

# %%
using StatsBase
g!(X, p) = (A[sample(keys(A), round(Int, p*length(A)); replace=false)] .= 0; A)
g!(copy(A), 0.25)

# %%
