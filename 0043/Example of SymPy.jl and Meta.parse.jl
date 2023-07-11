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
#     display_name: Julia 1.9.2
#     language: julia
#     name: julia-1.9
# ---

# %%
using SymPy
using QuadGK
using BenchmarkTools

@vars x real=true
@vars a positive=true

f(a, x) = √(a^2 + x^2)
expr_sym = sympy.Integral(f(a, x), x).doit()

# %%
str_sym = string(expr_sym)

# %%
expr_julia = Meta.parse("F(a, x) = $expr_sym")

# %%
@eval $expr_julia

# %%
f(x) = f(2, x)
F_quadgk(t) = quadgk(f, 0, t)[1]
F(t) = F(2, t) - F(2, 0)
t = 1:10^4
A = @btime F_quadgk.(t)
B = @btime F.(t);
A ≈ B

# %%
