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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
function repeatlength(A, i)
    k = 1
    @inbounds while i + k ≤ lastindex(A) && isequal(A[i+k], A[i])
        k += 1
    end
    k
end

"""
    RepVal(A)
Iterator generating the all (repeat length, value)-tuples in `A`.
"""
struct RepVal{T} A::T end
Base.IteratorSize(::Type{RepVal{T}}) where T = Base.SizeUnknown()
Base.eltype(x::RepVal) = Tuple{Int, eltype(x.A)} # (rep.len., value)

Base.iterate(x::RepVal) = iterate(x, firstindex(x.A))
function Base.iterate(x::RepVal, i::Int)
    i > lastindex(x.A) && return nothing
    k = repeatlength(x.A, i)
    (k, x.A[i]), i + k
end

maxrep_maxval(A) = maximum(RepVal(A))
negrep((k, v)) = (-k, v)
maxrep_minval(A) = negrep(minimum(negrep, RepVal(A)))

@doc RepVal

# %%
A = [1, 2, 3, 3, 3, 1, 1, 1, NaN, NaN, NaN, 2, 2, 3]
@show A
RepVal(A) |> collect

# %%
maxrep_maxval(A)

# %%
maxrep_minval(A)

# %%
@code_warntype iterate(RepVal(A), 1)

# %% [markdown]
# https://discourse.julialang.org/t/minimum-mode-problem-minimum-number-maximum-repetition/66749/2

# %%
function findn(x,y=zeros(Int,maximum(x)))
    m = 0
    fill!(y,0)
    @inbounds for i in x
        y[i] += 1
        m = max(y[i],m)
    end
    return findfirst(isequal(m),y)
end

# %%

# %%
using StatsBase, DataStructures

function findn_jling(x)
    cm = sort!(OrderedDict(countmap(x)); byvalue = true)
    last(cm.keys), last(cm.vals)
end

# %%

# %%
function minMode(c)   
    minVal = -1
    i = 0
    
    i = i+1
    for row in eachrow(c)
    row = collect(row)
    temp = filter(x -> x!=0, row)
    count = counter(temp)
    sortedCollection = sort(collect(count), by=x->x[2], rev=true)
    minVal = sortedCollection[1][1]
    repetitions = sortedCollection[1][2]
    for (key, value) in sortedCollection
        if(value == repetitions)
            if(key < minVal)
                minVal = key
            end
        end 
        end
    end
    
    return minVal
end

# %%
using StatsBase
swap_negval((val, rep)) = (rep, -val)
inv_swap_negval((rep, val)) = (-val, rep)
minmode(X) = inv_swap_negval(maximum(swap_negval, countmap(X)))k

# %%
X = rand(1:10000:10^8, 10^6)
X'

# %%
using BenchmarkTools

# %%
@btime maxrep_minval($X)

# %%
@btime findn($X)

# %%
@btime minMode($X)

# %%
@btime minmode($X)

# %%
@btime findn_jling($X)

# %%
