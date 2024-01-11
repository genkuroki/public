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

# %% [markdown]
# https://x.com/ziku1kokyu/status/1745293925014003788

# %%
using LinearAlgebra
dot2(x) = dot(x, x)

function answer(x::AbstractVector)
    n = length(x)
    c = 0
    for i in 1:n, j in i+1:n
        d2 = dot2(x[i] - x[j])
        c += d2 == isqrt(d2)^2
    end
    c
end

function answer(io::IO)
    n, d = parse.(Int, split(readline(io)))
    x = [parse.(Int, split(readline(io))) for _ in 1:n]
    c = answer(x)
    println(c)
end

# %%
problem = """
3 4
-3 7 8 2
-12 1 10 2
-2 8 9 3
"""
answer(IOBuffer(problem))

# %%
problem = """
3 4
-3 7 8 2
-12 1 10 2
-2 8 9 3
"""
write("problem.txt", problem)
open(io -> answer(io), "problem.txt")

# %%
problem = """
5 1
1
2
3
4
5
"""
answer(IOBuffer(problem))

# %%
answer(stdin)

# %%
