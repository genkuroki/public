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
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# https://twitter.com/servantprime/status/1639387897790083075
#
# <img src="73B62EC9-41F9-4E30-B6FC-3664F7D461AE.jpeg" width=60%>

# %%
using Primes

n = 10^9
P = primes(n)
@show length(P)
@show prod((p^2 + 1) / (p^2 - 1) for p in P);

# %%
a = prod(1 - 1/p^2 for p in P)
b = prod(1 + 1/p^2 for p in P)
c = prod(1 - 1/p^4 for p in P)
@show a b a*b c 6/π^2 90/π^4 b/a (90/6)/6;

# %%
