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
using SymPy
@vars x y z

# %%
a = Sym(2)^(1/Sym(3))

# %%
b = Sym(3)^(1/Sym(5))

# %%
θ = a + b

# %%
sympy.minimal_polynomial(θ, x)

# %%
f = sympy.minimal_polynomial(a, x)

# %%
g = sympy.minimal_polynomial(b, x)

# %%
h = -f(z - x).expand()

# %%
h1 = g - (x^2 + 3z*x + 6z^2)*h |> expand

# %%
r = (10z^3 - 2)^2*h - ((10z^3 - 2)*x - 15z^4 + 12z)*h1 |> expand

# %%
Denom, Numer = sympy.Poly(r, x).coeffs() |> C -> [C[1], -C[2]]

# %%
Numer / Denom

# %%
denom = Denom(z => θ) |> expand |> simplify

# %%
numer = Numer(z => θ) |> expand |> simplify

# %%
c = numer/denom |> simplify

# %%
c / b |> simplify

# %%
