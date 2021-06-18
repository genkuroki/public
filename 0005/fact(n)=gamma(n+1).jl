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
using Memoization
using SpecialFunctions
fact(n) = gamma(n + 1)

# %%
fact(100)

# %%
factorial(big(100)) |> Float64

# %%
coeff_sin(n) = Tuple((-1)^k/fact(2k+1) for k in 0:n-1)
taylor_odd(x, c) = x*evalpoly(x^2, c)

# %%
c = coeff_sin(8)
x = π/4
@show sin(x)
[(n, taylor_odd(x, coeff_sin(n))) for n in 1:9]

# %%
