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
#     display_name: Julia 1.7.0-beta3
#     language: julia
#     name: julia-1.7
# ---

# %%
using BenchmarkTools

f(n) = 1 + sum(abs2, 1:n)
g(n) = 1 + sum(abs2(k) for k in 1:n)

@show VERSION
@show f(10^6) == g(10^6)
print("f(10^6):")
@btime f(10^6)
print("g(10^6):")
@btime g(10^6)

# %%
using BenchmarkTools

f(n) = 1 + sum(abs2, 1:n)
g(n) = 1 + sum(abs2(k) for k in 1:n)

@show VERSION
@show f(10^6) == g(10^6)
print("f(10^6):")
@btime f(10^6)
print("g(10^6):")
@btime g(10^6)

# %%
import Base.Broadcast.preprocess_args
import Base.Broadcast.preprocess
@inline preprocess_args(dest, args::Tuple) = (Base.Broadcast.preprocess(dest, args[1]), Base.Broadcast.preprocess_args(dest, Base.tail(args))...)
@inline preprocess_args(dest, args::Tuple{Any}) = (Base.Broadcast.preprocess(dest, args[1]),)
@inline preprocess_args(dest, args::Tuple{}) = ()

using BenchmarkTools

f(n) = 1 + sum(abs2, 1:n)
g(n) = 1 + sum(abs2(k) for k in 1:n)

@show VERSION
@show f(10^6) == g(10^6)
print("f(10^6):")
@btime f(10^6)
print("g(10^6):")
@btime g(10^6)

# %%
