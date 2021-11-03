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
using BenchmarkTools
using CUDA
using Random

# %%
mcpi_cu_sum(n, t::Type{T}=Float32) where T = 4count(≤(1), sum(x->x^2, CUDA.rand(T, 2, n); dims=1))/n

@time mcpi_cu_sum(10^8)
@time mcpi_cu_sum(10^8)
@time mcpi_cu_sum(10^8)

# %%
@time mcpi_cu_sum(10^8, Float64)
@time mcpi_cu_sum(10^8, Float64)
@time mcpi_cu_sum(10^8, Float64)

# %%
mcpi_cu_count(n, t::Type{T}=Float32) where T = 4count(x -> x^2 + rand(T)^2 ≤ 1, CUDA.rand(T, n))/n

@time mcpi_cu_count(10^8)
@time mcpi_cu_count(10^8)
@time mcpi_cu_count(10^8)

# %%
mcpi_cu_count(n, t::Type{T}=Float32) where T = 4count(x -> x^2 + rand(T)^2 ≤ 1, CUDA.rand(T, n))/n

@time mcpi_cu_count(10^8, Float64)
@time mcpi_cu_count(10^8, Float64)
@time mcpi_cu_count(10^8, Float64)

# %%
function mcpi(n, t::Type{T}=Float32) where T
    rng = MersenneTwister()
    4count(_ -> rand(rng, T)^2 + rand(rng, T)^2 ≤ 1, 1:n)/n
end

@time mcpi(10^8)
@time mcpi(10^8)
@time mcpi(10^8)

# %%
@time mcpi(10^8, Float64)
@time mcpi(10^8, Float64)
@time mcpi(10^8, Float64)

# %%
@btime mcpi(10^8)
@btime mcpi_cu_sum(10^8)
@btime mcpi_cu_count(10^8)

# %%
@btime mcpi(10^8, Float64)
@btime mcpi_cu_sum(10^8, Float64)
@btime mcpi_cu_count(10^8, Float64)

# %%
