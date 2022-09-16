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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# ガンマ函数は次のように定義される:
#
# $$
# \Gamma(s) = \int_0^\infty e^{-x} x^{s-1}\,dx
# $$

# %%
using SpecialFunctions
using Plots
plot(loggamma, 0.01, 10.0; legend=:topleft, label="loggamma") 

# %%
plot(gamma, 0.5, 3.0; legend=:topleft, label="gamma")

# %%
