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
function mcmcpi(L=10^9)
    c = 0
    for i in 1:L
        c += rand()^2 + rand()^2 ≤ 1
    end
    4c/L
end

@time mcmcpi()
@time mcmcpi()
@time mcmcpi()

# %%
function mcmcpi_threads(L=10^9)
    c = zeros(Int, Threads.nthreads())
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        c[tid] += rand()^2 + rand()^2 ≤ 1
    end
    4sum(c)/L
end

@time mcmcpi_threads()
@time mcmcpi_threads()
@time mcmcpi_threads()

# %%
using LoopVectorization

# The result shall be 4.0 or 0.0.
function mcmcpi_turbo_incorrect(L=10^9)
    c = 0
    @turbo for i in 1:L
        c += rand()^2 + rand()^2 ≤ 1
    end
    4c/L
end

@time mcmcpi_turbo_incorrect()
@time mcmcpi_turbo_incorrect()
@time mcmcpi_turbo_incorrect()

# %%
using LoopVectorization

rand_is_in_unit_circle(i) = rand()^2 + rand()^2 ≤ 1

function mcmcpi_turbo(L=10^9)
    c = 0
    @turbo for i in 1:L
        c += rand_is_in_unit_circle(i)
    end
    4c/L
end

@time mcmcpi_turbo()
@time mcmcpi_turbo()
@time mcmcpi_turbo()

# %%
