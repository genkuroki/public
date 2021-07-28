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
using SymPy
@syms z
@syms x::real
@syms a::positive

# %%
expr = cos(z)^2 + sin(z)^2

# %%
expr.simplify()

# %%
Eq(expr, expr.simplify())

# %%
G = sympy.Integral(exp(-x^2/a), (x, -oo, oo))

# %%
Eq(G, G.doit())

# %%
