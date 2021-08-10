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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using BenchmarkTools
b = @benchmark sum(rand()^2 + rand()^2 ≤ 1 for _ in 1:10^6)/10^6

# %%
dump(b)

# %%
methodswith(BenchmarkTools.Trial)

# %%
b.times

# %%
