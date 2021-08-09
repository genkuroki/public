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
using Test
L = 20;

# %%
function f(L)
    a = zeros(Int, 2^L)
    for j in 0:2^L-1
        for i in 0:L-1
            a[begin + (j ⊻ (1 << i))] += 1
        end
    end
    a
end
@test f(L) == fill(L, 2^L)

# %%
function f_ill_multithreaded(L)
    a = zeros(Int, 2^L)
    Threads.@threads for j in 0:2^L-1
        for i in 0:L-1
            a[begin + (j ⊻ (1 << i))] += 1
        end
    end
    a
end
@test f_ill_multithreaded(L) == fill(L, 2^L)

# %%
function f_well_multithreaded(L)
    a = zeros(Int, 2^L)
    for i in 0:L-1
        Threads.@threads for j in 0:2^L-1
            a[begin + (j ⊻ (1 << i))] += 1
        end
    end
    a
end
@test f_well_multithreaded(L) == fill(L, 2^L)

# %%
