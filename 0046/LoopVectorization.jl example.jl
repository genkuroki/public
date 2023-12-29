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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
function mcpi(L=10^9)
    c = 0
    for i in 1:L
        c += rand()^2 + rand()^2 ≤ 1
    end
    4c/L
end

@time mcpi()
@time mcpi()
@time mcpi()

# %%
function mcpi_threads(L=10^9)
    c = zeros(Int, Threads.nthreads())
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        c[tid] += rand()^2 + rand()^2 ≤ 1
    end
    4sum(c)/L
end

@time mcpi_threads()
@time mcpi_threads()
@time mcpi_threads()

# %%
using LoopVectorization

# The result shall be 4.0 or 0.0.
function mcpi_turbo_incorrect(L=10^9)
    c = 0
    @turbo for i in 1:L
        c += rand()^2 + rand()^2 ≤ 1
    end
    4c/L
end

@time mcpi_turbo_incorrect()
@time mcpi_turbo_incorrect()
@time mcpi_turbo_incorrect()

# %%
using LoopVectorization

rand_is_in_unit_circle(i) = rand()^2 + rand()^2 ≤ 1

function mcpi_turbo(L=10^9)
    c = 0
    @turbo for i in 1:L
        c += rand_is_in_unit_circle(i)
    end
    4c/L
end

@time mcpi_turbo()
@time mcpi_turbo()
@time mcpi_turbo()

# %%
using CUDA
using Statistics

square(x) = x^2

function mcpi_cuda_count(L=10^8, t::Type{T}=Float32) where T
    4count(≤(1), sum(square, CUDA.rand(T, 2, L); dims=1)) / L
end

@time mean(mcpi_cuda_count() for _ in 1:10)
@time mean(mcpi_cuda_count() for _ in 1:10)
@time mean(mcpi_cuda_count() for _ in 1:10)

# %%
