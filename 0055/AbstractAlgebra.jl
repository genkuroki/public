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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# https://colab.research.google.com/github/genkuroki/public/blob/main/0055/AbstractAlgebra.ipynb

# %%
haskey(ENV, "COLAB_GPU") && (import Pkg; Pkg.add("AbstractAlgebra"))
using AbstractAlgebra

# %%
P, a = polynomial_ring(QQ, :a)

# %%
R, = residue_ring(P, a^2 - 11)

# %%
(4a)^2 |> R # √176 = 4a

# %%
f(x) = x^4 + 4a*x + 2 # = x^4 + √176 x + 2
f(2a) |> R

# %%
