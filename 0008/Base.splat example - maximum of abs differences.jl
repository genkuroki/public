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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://discourse.julialang.org/t/is-there-a-maximum-f-op-itrs-in-julia/63868/3

# %%
using BenchmarkTools

a = randn(1000, 1000)
b = randn(1000, 1000)

f(a, b) = maximum(abs, x - y for (x, y) in zip(a, b))
g(a, b) = maximum(Base.splat(abs∘-), zip(a, b))
h(a, b) = mapreduce(abs∘-, max, a, b)

@show f(a, b) == g(a, b) == h(a, b)
@btime f($a, $b)
@btime g($a, $b)
@btime h($a, $b)

# %%
