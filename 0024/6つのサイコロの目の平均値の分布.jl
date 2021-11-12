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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using Statistics

# %%
X = [mean(rand(1:6, 6)) for _ in 1:10^6]
mean(==(3.5), X)

# %%
mean(≈(3.5), X)

# %%
X = [sum(X) for X in Iterators.product(fill(1:6, 6)...)]
count(==(21), X)//6^6

# %%
mean(==(21), X)

# %%
