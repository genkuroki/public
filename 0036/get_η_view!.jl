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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %%
using BenchmarkTools
macro atime(expr) :(@btime $expr samples=1 evals=1) end

# %%
function get_η_(T)
    D = ndims(T)
    η = reverse(T)
    for d in 1:D
        cumsum!(η, η, dims=d)
    end
    reverse!(η)
end

# %%
function get_η_view!(T)
    D = ndims(T)
    η = view(T, reverse.(axes(T))...)
    for d in 1:D
        cumsum!(η, η, dims=d)
    end
    parent(η)
end

get_η_view(T) = get_η_view!(copy(T))

# %%
T = rand(1:8, 5, 4, 3)
A = get_η_(T)
B = get_η_view(T)
A == B

# %%
T = @atime rand(1:8, 400, 400, 400);

# %%
A = @atime get_η_($T);

# %%
B = @atime get_η_view($T);

# %%
A == B

# %%
tmp = @time copy(T)
C = @time get_η_view!(tmp)
A == C

# %%
