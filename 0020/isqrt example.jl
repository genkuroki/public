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
ENV["LINES"] = 1000

# %%
n = 100
[(a, b, isqrt(a^2+b^2)) for a in 3:n for b in a+1:a^2÷2
    if gcd(a, b) == 1 && a^2 + b^2 == isqrt(a^2 + b^2)^2]

# %%
n = 100
[(a, b, isqrt(a^2+b^2)) for a in 3:n for b in a+1:max(n, a^2÷2)
    if gcd(a, b) == 1 && a^2 + b^2 == isqrt(a^2 + b^2)^2 && b != a^2÷2]

# %%
