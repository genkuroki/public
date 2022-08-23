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
A = rand(0:9, 2, 3, 4, 5)
perm = (4, 3, 1, 2)
@code_warntype PermutedDimsArray(A, perm)

# %%
@which PermutedDimsArray(zeros(1, 2, 3, 4), (4, 2, 1, 3))

# %%
module O

struct PermutedDimsArray{T,N,AA<:AbstractArray{T,N},NT<:Dims{N}} <: AbstractArray{T,N}
    parent::AA
    perm::NT
    iperm::NT
end

function PermutedDimsArray(data::AbstractArray{T,N}, perm::Dims{N}) where {T,N}
    iperm = invperm(perm)
    PermutedDimsArray{T,N,typeof(data),typeof(perm)}(data, perm, iperm)
end
function PermutedDimsArray(data::AbstractArray{T,N}, perm) where {T,N}
    PermutedDimsArray(data, Dims{N}(perm))
end

Base.parent(A::PermutedDimsArray) = A.parent
Base.size(A::PermutedDimsArray{T,N}) where {T,N} = genperm(size(parent(A)), A.perm)
Base.axes(A::PermutedDimsArray{T,N}) where {T,N} = genperm(axes(parent(A)), A.perm)

Base.similar(A::PermutedDimsArray, T::Type, dims::Base.Dims) = similar(parent(A), T, dims)

Base.unsafe_convert(::Type{Ptr{T}}, A::PermutedDimsArray{T}) where {T} = Base.unsafe_convert(Ptr{T}, parent(A))

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

end

# %%
A = rand(0:9, 2, 3, 4, 5)
perm = (4, 3, 1, 2)
@code_warntype O.PermutedDimsArray(A, perm)

# %%
A = rand(2, 3, 4, 5)
perm = (4, 3, 1, 2)
B = @atime PermutedDimsArray($A, $perm)
C = @atime O.PermutedDimsArray($A, $perm)
B == C

# %%
A = rand(2, 3, 4, 5)
perm = [4, 3, 1, 2]
B = @atime PermutedDimsArray($A, $perm)
C = @atime O.PermutedDimsArray($A, $perm)
B == C

# %%
