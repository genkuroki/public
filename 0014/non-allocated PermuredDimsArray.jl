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

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/073af4a7415c6e9f15df29b6cf3732050f4ce7c8/base/permuteddimsarray.jl

# %%
# This file is a part of Julia. License is MIT: https://julialang.org/license

module My

import Base: permutedims, permutedims!
# export PermutedDimsArray

# Some day we will want storage-order-aware iteration, so put perm in the parameters
struct PermutedDimsArray{T,N,AA<:AbstractArray{T,N}} <: AbstractArray{T,N}
    parent::AA
    perm::NTuple{N,Int}
    iperm::NTuple{N,Int}

    function PermutedDimsArray{T,N,AA}(data::AA, perm::NTuple{N,Int}, iperm::NTuple{N,Int}) where {T,N,AA<:AbstractArray{T,N}}
        isperm(perm) || throw(ArgumentError(string(perm, " is not a valid permutation of dimensions 1:", N)))
        all(iperm[perm[d]]==d for d in 1:N) || throw(ArgumentError(string(perm, " and ", iperm, " must be inverses")))
        new(data, perm, iperm)
    end
end

"""
    PermutedDimsArray(A, perm) -> B
Given an AbstractArray `A`, create a view `B` such that the
dimensions appear to be permuted. Similar to `permutedims`, except
that no copying occurs (`B` shares storage with `A`).
See also [`permutedims`](@ref), [`invperm`](@ref).
# Examples
```jldoctest
julia> A = rand(3,5,4);
julia> B = PermutedDimsArray(A, (3,1,2));
julia> size(B)
(4, 3, 5)
julia> B[3,1,2] == A[1,2,3]
true
```
"""
function PermutedDimsArray(data::AbstractArray{T,N}, perm) where {T,N}
    length(perm) == N || throw(ArgumentError(string(perm, " is not a valid permutation of dimensions 1:", N)))
    iperm = invperm(perm)
    PermutedDimsArray{T,N,typeof(data)}(data, (perm...,), (iperm...,))
end

Base.parent(A::PermutedDimsArray) = A.parent
Base.size(A::PermutedDimsArray) = genperm(size(parent(A)), A.perm)
Base.axes(A::PermutedDimsArray) = genperm(axes(parent(A)), A.perm)

Base.similar(A::PermutedDimsArray, T::Type, dims::Base.Dims) = similar(parent(A), T, dims)

Base.unsafe_convert(::Type{Ptr{T}}, A::PermutedDimsArray{T}) where {T} = Base.unsafe_convert(Ptr{T}, parent(A))

# It's OK to return a pointer to the first element, and indeed quite
# useful for wrapping C routines that require a different storage
# order than used by Julia. But for an array with unconventional
# storage order, a linear offset is ambiguous---is it a memory offset
# or a linear index?
Base.pointer(A::PermutedDimsArray, i::Integer) = throw(ArgumentError("pointer(A, i) is deliberately unsupported for PermutedDimsArray"))

function Base.strides(A::PermutedDimsArray{T,N}) where {T,N}
    s = strides(parent(A))
    ntuple(d->s[A.perm[d]], Val(N))
end
Base.elsize(::Type{<:PermutedDimsArray{<:Any, <:Any, P}}) where {P} = Base.elsize(P)

@inline function Base.getindex(A::PermutedDimsArray{T,N}, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds val = getindex(A.parent, genperm(I, A.iperm)...)
    val
end
@inline function Base.setindex!(A::PermutedDimsArray{T,N}, val, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds setindex!(A.parent, val, genperm(I, A.iperm)...)
    val
end

@inline genperm(I::NTuple{N,Any}, perm::Dims{N}) where {N} = ntuple(d -> I[perm[d]], Val(N))
@inline genperm(I, perm::AbstractVector{Int}) = genperm(I, (perm...,))

"""
    permutedims(A::AbstractArray, perm)
Permute the dimensions of array `A`. `perm` is a vector or a tuple of length `ndims(A)`
specifying the permutation.
See also [`permutedims!`](@ref), [`PermutedDimsArray`](@ref), [`transpose`](@ref), [`invperm`](@ref).
# Examples
```jldoctest
julia> A = reshape(Vector(1:8), (2,2,2))
2×2×2 Array{Int64, 3}:
[:, :, 1] =
 1  3
 2  4
[:, :, 2] =
 5  7
 6  8
julia> permutedims(A, (3, 2, 1))
2×2×2 Array{Int64, 3}:
[:, :, 1] =
 1  3
 5  7
[:, :, 2] =
 2  4
 6  8
julia> B = randn(5, 7, 11, 13);
julia> perm = [4,1,3,2];
julia> size(permutedims(B, perm))
(13, 5, 11, 7)
julia> size(B)[perm] == ans
true
```
"""
function permutedims(A::AbstractArray, perm)
    dest = similar(A, genperm(axes(A), perm))
    permutedims!(dest, A, perm)
end

"""
    permutedims(m::AbstractMatrix)
Permute the dimensions of the matrix `m`, by flipping the elements across the diagonal of
the matrix. Differs from `LinearAlgebra`'s [`transpose`](@ref) in that the
operation is not recursive.
# Examples
```jldoctest; setup = :(using LinearAlgebra)
julia> a = [1 2; 3 4];
julia> b = [5 6; 7 8];
julia> c = [9 10; 11 12];
julia> d = [13 14; 15 16];
julia> X = [[a] [b]; [c] [d]]
2×2 Matrix{Matrix{Int64}}:
 [1 2; 3 4]     [5 6; 7 8]
 [9 10; 11 12]  [13 14; 15 16]
julia> permutedims(X)
2×2 Matrix{Matrix{Int64}}:
 [1 2; 3 4]  [9 10; 11 12]
 [5 6; 7 8]  [13 14; 15 16]
julia> transpose(X)
2×2 transpose(::Matrix{Matrix{Int64}}) with eltype Transpose{Int64, Matrix{Int64}}:
 [1 3; 2 4]  [9 11; 10 12]
 [5 7; 6 8]  [13 15; 14 16]
```
"""
permutedims(A::AbstractMatrix) = permutedims(A, (2,1))

"""
    permutedims(v::AbstractVector)
Reshape vector `v` into a `1 × length(v)` row matrix.
Differs from `LinearAlgebra`'s [`transpose`](@ref) in that
the operation is not recursive.
# Examples
```jldoctest; setup = :(using LinearAlgebra)
julia> permutedims([1, 2, 3, 4])
1×4 Matrix{Int64}:
 1  2  3  4
julia> V = [[[1 2; 3 4]]; [[5 6; 7 8]]]
2-element Vector{Matrix{Int64}}:
 [1 2; 3 4]
 [5 6; 7 8]
julia> permutedims(V)
1×2 Matrix{Matrix{Int64}}:
 [1 2; 3 4]  [5 6; 7 8]
julia> transpose(V)
1×2 transpose(::Vector{Matrix{Int64}}) with eltype Transpose{Int64, Matrix{Int64}}:
 [1 3; 2 4]  [5 7; 6 8]
```
"""
permutedims(v::AbstractVector) = reshape(v, (1, length(v)))

"""
    permutedims!(dest, src, perm)
Permute the dimensions of array `src` and store the result in the array `dest`. `perm` is a
vector specifying a permutation of length `ndims(src)`. The preallocated array `dest` should
have `size(dest) == size(src)[perm]` and is completely overwritten. No in-place permutation
is supported and unexpected results will happen if `src` and `dest` have overlapping memory
regions.
See also [`permutedims`](@ref).
"""
function permutedims!(dest, src::AbstractArray, perm)
    Base.checkdims_perm(dest, src, perm)
    P = PermutedDimsArray(dest, invperm(perm))
    _copy!(P, src)
    return dest
end

function Base.copyto!(dest::PermutedDimsArray{T,N}, src::AbstractArray{T,N}) where {T,N}
    checkbounds(dest, axes(src)...)
    _copy!(dest, src)
end
Base.copyto!(dest::PermutedDimsArray, src::AbstractArray) = _copy!(dest, src)

function _copy!(P::PermutedDimsArray, src)
    # If dest/src are "close to dense," then it pays to be cache-friendly.
    # Determine the first permuted dimension
    d = 0  # d+1 will hold the first permuted dimension of src
    while d < ndims(src) && P.perm[d+1] == d+1
        d += 1
    end
    if d == ndims(src)
        copyto!(parent(P), src) # it's not permuted
    else
        R1 = CartesianIndices(axes(src)[1:d])
        d1 = findfirst(isequal(d+1), perm)::Int  # first permuted dim of dest
        R2 = CartesianIndices(axes(src)[d+2:d1-1])
        R3 = CartesianIndices(axes(src)[d1+1:end])
        _permutedims!(P, src, R1, R2, R3, d+1, d1)
    end
    return P
end

@noinline function _permutedims!(P::PermutedDimsArray, src, R1::CartesianIndices{0}, R2, R3, ds, dp)
    ip, is = axes(src, dp), axes(src, ds)
    for jo in first(ip):8:last(ip), io in first(is):8:last(is)
        for I3 in R3, I2 in R2
            for j in jo:min(jo+7, last(ip))
                for i in io:min(io+7, last(is))
                    @inbounds P[i, I2, j, I3] = src[i, I2, j, I3]
                end
            end
        end
    end
    P
end

@noinline function _permutedims!(P::PermutedDimsArray, src, R1, R2, R3, ds, dp)
    ip, is = axes(src, dp), axes(src, ds)
    for jo in first(ip):8:last(ip), io in first(is):8:last(is)
        for I3 in R3, I2 in R2
            for j in jo:min(jo+7, last(ip))
                for i in io:min(io+7, last(is))
                    for I1 in R1
                        @inbounds P[I1, i, I2, j, I3] = src[I1, i, I2, j, I3]
                    end
                end
            end
        end
    end
    P
end

function Base._mapreduce_dim(f, op, init::Base._InitialValue, A::PermutedDimsArray, dims::Colon)
    Base._mapreduce_dim(f, op, init, parent(A), dims)
end

function Base.mapreducedim!(f, op, B::AbstractArray{T,N}, A::PermutedDimsArray{T,N}) where {T,N}
    C = PermutedDimsArray{T,N,typeof(B)}(B, B.perm, B.iperm) # make the inverse permutation for the output
    Base.mapreducedim!(f, op, C, parent(A))
    B
end

function Base.showarg(io::IO, A::PermutedDimsArray{T,N}, toplevel) where {T,N}
    print(io, "PermutedDimsArray(")
    Base.showarg(io, parent(A), false)
    print(io, ", ", A.perm, ')')
    toplevel && print(io, " with eltype ", eltype(A))
    return nothing
end

end

# %%
using BenchmarkTools
x = 100(1:4) .+ 10(1:3)' .+ reshape(1:2, 1, 1, :)
@btime p = My.PermutedDimsArray($x, (2, 3, 1))

# %%
@btime q = PermutedDimsArray($x, (2, 3, 1))

# %%
