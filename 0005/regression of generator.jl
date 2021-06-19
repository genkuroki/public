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

f(n) = sum(abs2, 1:n)
g(n) = sum(abs2(x) for x in 1:n)

@show VERSION
@btime f(10^6)
@btime g(10^6)

# %%
using BenchmarkTools

f(n) = sum(abs2, 1:n)
g(n) = sum(abs2(x) for x in 1:n)

@show VERSION
@btime f(10^6)
@btime g(10^6)

# %%
using BenchmarkTools

f(n) = sum(abs2, 1:n)
g(n) = sum(abs2(x) for x in 1:n)

@show VERSION
@btime f(10^6)
@btime g(10^6)

# %%
