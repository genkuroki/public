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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using SparseArrays

nzptr(A::SparseArrays.AbstractSparseMatrixCSC, col::Integer) = view(rowvals(A), nzrange(A, col))

using Random, LinearAlgebra
Random.seed!(12345678)

n = 20000
p = 2/n
@time A = Symmetric(sprand(Bool, n, n, p)) # test adjacency matrix
@time B = sparse(A)
@time C = nzptr.(Ref(B), axes(B, 2)) # adjacency list

# %%
B[:, 4]

# %%
Random.seed!(12345678)

n = 20000
p = 0.1
@time A = Symmetric(sprand(Bool, n, n, p)) # test adjacency matrix
@time B = sparse(A)
@time C = nzptr.(Ref(B), axes(B, 2)) # adjacency list

# %%
B[:, 1]

# %%
