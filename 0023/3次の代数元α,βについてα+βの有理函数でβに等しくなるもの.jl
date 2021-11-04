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

# %% tags=[]
using SymPy
@vars X t a b c p q r α β

# %% tags=[]
f = X^3 + a*X^2 + b*X + c

# %%
g = X^3 + p*X^2 + q*X + r

# %%
h = -f(X => t - X).expand()

# %%
h1 = g - h

# %%
h2 = (a + p + 3t)^2*h - ((a + p + 3t)*X + (-a^2 - a*p - 4a*t + b - 3p*t - q - 6t^2))*h1 |> expand |> simplify

# %%
Denom, Numer = sympy.Poly(h2, X).coeffs() |> C -> [C[1], -C[2]];

# %%
DENOM = sympy.Poly(Denom, t)

# %%
NUMER = sympy.Poly(Numer, t)

# %%
βdenom = (β * Denom(t => α + β)).expand()(α^5 => α^5 - α^2*f(X=>α), β^5 => β^5 - β^2*g(X=>β)).expand()(α^4 => α^4 - α*f(X=>α), β^4 => β^4 - β*g(X=>β)).expand().(α^3 => α^3 - f(X=>α), β^3 => β^3 - g(X=>β)).expand()

# %%
numer = Numer(t => α + β).expand()(α^5 => α^5 - α^2*f(X=>α), β^5 => β^5 - β^2*g(X=>β)).expand()(α^4 => α^4 - α*f(X=>α), β^4 => β^4 - β*g(X=>β)).expand().(α^3 => α^3 - f(X=>α), β^3 => β^3 - g(X=>β)).expand()

# %%
βdenom == numer

# %%
