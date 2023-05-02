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

# %%
using Primes

# %% tags=[]
N = 20000
A = [(p, p+2, p+6, p+8) for p in primes(6, N) if isprime(p+2) && isprime(p+6) && isprime(p+8)]

# %%
[mod(sum(a), 60) == 0 for a in A] |> all

# %% [markdown]
# 問題: $p$, $p+2$, $p+6$, $p+8$ が7以上の素数ならば, それら4つの素数の和が60で割り切れることを示せ.

# %%
