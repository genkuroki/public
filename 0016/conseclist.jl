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
# Subset of https://gist.github.com/genkuroki/8c7e9e484a5877d9747a4819bbab645d

# %%
function conseclist_forloop(x)
    b, e = firstindex(x), lastindex(x)
    s = b - 1
    S = [s]
    while s < e
        t = s + 1
        @inbounds while t < e
            x[t+1] != x[t] + 1 && break 
            t += 1
        end
        push!(S, t)
        s = t
    end
    [@view(x[S[i]+1:S[i+1]]) for i in 1:length(S)-1]
end

# %%
function conseclist_bitvector(x)
    b, e = firstindex(x), lastindex(x)
    @views idxs = x[b+1:e] .- x[b:e-1] .!= 1
    S = [b-1; (b:e-1)[idxs]; e]
    [@view(x[S[i]+1:S[i+1]]) for i in 1:length(S)-1]
end

# %%
module LazyConsecLists

struct ConsecList{T<:AbstractVector{<:Integer}} <:
        AbstractVector{SubArray{eltype(T), 1, T, Tuple{UnitRange{Int64}}, true}}
    a::T
    S::Vector{Int}
end
function ConsecList(x)
    b, e = firstindex(x), lastindex(x)
    idx = @views (x[b+1:e] .- x[b:e-1]) .!= 1
    S = [b-1; (b:e-1)[idx]; e]
    ConsecList{typeof(x)}(x, S)
end

Base.length(x::ConsecList) = length(x.S) - 1
Base.size(x::ConsecList) = (length(x),)
function Base.eltype(x::ConsecList{T}) where T<:AbstractVector{<:Integer}
    SubArray{eltype(T), 1, T, Tuple{UnitRange{Int64}}, true}
end

Base.getindex(x::ConsecList, i::Integer) = @view(x.a[x.S[i]+1:x.S[i+1]])
Base.getindex(x::ConsecList, r::AbstractRange) = [x[i] for i in r]
Base.getindex(x::ConsecList, ::Colon) = collect(x)

end

# %%
using BenchmarkTools
A = unique(sort(rand(1:10^4, 10^4)))

# %%
@show conseclist_forloop(A) == conseclist_bitvector(A) == LazyConsecLists.ConsecList(A)
@show length(conseclist_forloop(A)) == length(conseclist_bitvector(A)) == length(LazyConsecLists.ConsecList(A))
@btime conseclist_forloop($A)
@btime conseclist_bitvector($A)
@btime collect(LazyConsecLists.ConsecList($A))

# %%
@btime sum(length, conseclist_forloop($A))
@btime sum(length, conseclist_bitvector($A))
@btime sum(length, LazyConsecLists.ConsecList($A))

# %%
