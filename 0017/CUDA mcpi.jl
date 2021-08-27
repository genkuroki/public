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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
using BenchmarkTools
using Random
function mcpi_f32(n)
    rng = MersenneTwister()
    4count(_ -> rand(rng, Float32)^2 + rand(rng, Float32)^2 ≤ 1, 1:n)/n
end
@show mcpi_f32(10^8)
@btime mcpi_f32(10^8);

# %%
using BenchmarkTools
using CUDA
mcpi_f32_cu(n) = 4count(x -> x^2 + rand(Float32)^2 ≤ 1, CUDA.rand(n))/n
@show mcpi_f32_cu(10^8)
@btime mcpi_f32_cu(10^8);

# %%
using BenchmarkTools
using Random
function mcpi_f64(n)
    rng = MersenneTwister()
    4count(_ -> rand(rng)^2 + rand(rng)^2 ≤ 1, 1:n)/n
end
@show mcpi_f64(10^8)
@btime mcpi_f64(10^8);

# %%
using BenchmarkTools
using CUDA
mcpi_f64_cu(n) = 4count(x -> x^2 + rand(Float64)^2 ≤ 1, CUDA.rand(Float64, n))/n
@show mcpi_f64_cu(10^8)
@btime mcpi_f64_cu(10^8);

# %%
CUDA.versioninfo()

# %%
