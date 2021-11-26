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
#     display_name: Julia 1.6.4
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://twitter.com/genkuroki/status/1463113855299584007

# %%
using SymPy
@vars x w μ σ² f real=true
N(x, μ, σ²) = 1/√(2oftype(x, π)*σ²) * exp(-(x - μ)^2/(2σ²))
N(x, μ, σ²)

# %%
var"P(w)" = N(w, 0, 1)

# %%
var"P(x|w)" = N(x, w, 1)

# %%
var"P(x, w)" = var"P(x|w)" * var"P(w)" |> expand |> simplify

# %%
var"P(x)" = integrate(var"P(x, w)", (w, -oo, oo)) |> expand |> simplify

# %%
var"P(x)" == N(x, 0, 2)

# %%
var"P(w|x)" = var"P(x, w)"/var"P(x)" |> expand |> simplify

# %%
var"P(w|x)" == N(w, x/2, 1//2) |> expand |> simplify

# %%
var"∫L(f, w)P(w|x)dw" = integrate((f - w)^2 * var"P(w|x)", (w, -oo, oo)) |> expand |> simplify

# %%
var"∫L(f, w)P(w|x)dw" == (f - x/2)^2 + 1//2 |> expand

# %%
