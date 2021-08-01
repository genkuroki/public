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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
VERSION

# %%
indnearest_alloc(M, x) = argmin(@. abs(M - x)) # allocates!

# See https://github.com/JuliaLang/julia/blob/master/base/reduce.jl#L862
valindargmin(f, X) = mapfoldl(((i, x),) -> (f(x), i, x),
    ((v, i, x), (w, j, y)) -> v ≥ w ? (w, j, y) : (v, i, x), pairs(X))
valindargnearest(M, x) = valindargmin(m -> abs(m - x), M)
valindnearest(M, x) = ((val, ind, arg) = valindargnearest(M, x); (val, ind))
indargnearest(M, x) = ((val, ind, arg) = valindargnearest(M, x); (ind, arg))
indnearest(M, x) = ((val, ind, arg) = valindargnearest(M, x); ind)

# %%
f(x, y) = sin(x) * cos(x) + 0.01x * y
M = f.(2.5:0.2:3.7, (1:0.2:2)')

# %%
indnearest_alloc(M, 0.2)

# %%
valindargnearest(M, 0.2)

# %%
valindnearest(M, 0.2)

# %%
indargnearest(M, 0.2)

# %%
indnearest(M, 0.2)

# %%
using BenchmarkTools
@btime indnearest_alloc($M, 0.2)
@btime valindargnearest($M, 0.2)
@btime valindnearest($M, 0.2)
@btime indargnearest($M, 0.2)
@btime indnearest($M, 0.2);

# %%
