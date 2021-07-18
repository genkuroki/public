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

# %% [markdown]
# https://discourse.julialang.org/t/how-to-reduce-large-memory-allocation-by-lu-with-sparse-matrices/64791/7

# %%
VERSION

# %%
using BenchmarkTools
using Random; Random.seed!(4649373)
using LinearAlgebra
using SparseArrays
using SparseArrays: getcolptr
using SuiteSparse: decrement
using SuiteSparse.UMFPACK: UmfpackLU, UMFVTypes, UMFITypes, umfpack_numeric!

"""https://github.com/JuliaLang/SuiteSparse.jl/blob/master/src/umfpack.jl#L220-L278"""
function lu_org!(F::UmfpackLU, S::SparseMatrixCSC{<:UMFVTypes,<:UMFITypes}; check::Bool=true)
    zerobased = getcolptr(S)[1] == 0
    F.m = size(S, 1)
    F.n = size(S, 2)
    F.colptr = zerobased ? copy(getcolptr(S)) : decrement(getcolptr(S))
    F.rowval = zerobased ? copy(rowvals(S)) : decrement(rowvals(S))
    F.nzval = copy(nonzeros(S))

    umfpack_numeric!(F, reuse_numeric = false)
    check && (issuccess(F) || throw(LinearAlgebra.SingularException(0)))
    return F
end

"""A modification of https://github.com/JuliaLang/SuiteSparse.jl/blob/master/src/umfpack.jl#L220-L278"""
function lu_dot!(F::UmfpackLU, S::SparseMatrixCSC{<:UMFVTypes,<:UMFITypes}; check::Bool=true)
    zerobased = getcolptr(S)[1] == 0
    F.m = size(S, 1)
    F.n = size(S, 2)
    if zerobased
        F.colptr .= getcolptr(S)
        F.rowval .= rowvals(S)
    else
        F.colptr .= getcolptr(S) .- oneunit(eltype(S))
        F.rowval .= rowvals(S) .- oneunit(eltype(S))
    end
    F.nzval .= nonzeros(S)

    umfpack_numeric!(F, reuse_numeric = false)
    check && (issuccess(F) || throw(LinearAlgebra.SingularException(0)))
    return F
end

n = 2^10
X = I + sprandn(n, n, 1e-3)

F = lu(X)
G = lu(X)
@show lu(X).L == lu!(F, X).L
@show lu(X).U == lu!(F, X).U
@show lu!(F, X).L == lu_dot!(G, X).L
@show lu!(F, X).U == lu_dot!(G, X).U
println()
print("lu(X):        ")
@btime lu($X)
print("lu!(F, X):    ")
@btime lu!($F, $X)
print("lu_dot!(F, X):")
@btime lu_dot!($F, $X);

# %%
print("umfpack_numeric!(F, reuse_numeric = false):")
let F = F, S = X
    zerobased = getcolptr(S)[1] == 0
    F.m = size(S, 1)
    F.n = size(S, 2)
    if zerobased
        F.colptr .= getcolptr(S)
        F.rowval .= rowvals(S)
    else
        F.colptr .= getcolptr(S) .- oneunit(eltype(S))
        F.rowval .= rowvals(S) .- oneunit(eltype(S))
    end
    F.nzval .= nonzeros(S)

    @btime umfpack_numeric!(F, reuse_numeric = false)
end;

# %%
@which umfpack_numeric!(F, reuse_numeric = false)

# %% [markdown]
# https://github.com/JuliaLang/SuiteSparse.jl/blob/master/src/umfpack.jl#L384-L400

# %%
mutable struct Factor F::UmfpackLU end
factor = Factor(F)
factor.F === F

# %%
