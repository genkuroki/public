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
# * https://github.com/JuliaMath/Combinatorics.jl/blob/master/src/permutations.jl#L47
# * https://discourse.julialang.org/t/is-there-a-function-behaving-the-same-as-next-permutation-does-in-c/63451/13

# %%
module O

function next_perm!(v::AbstractVector)
    length(v) ≤ 1 && return false
    
    k = findlast(isless(v[i], v[i+1]) for i in firstindex(v):lastindex(v)-1)
    isnothing(k) && (reverse!(v); return false)
    i = k + firstindex(v)-1
    
    j = findlast(isless(v[i], v[j]) for j in i+1:lastindex(v)) + i
    v[i], v[j] = v[j], v[i]
    reverse!(v, i + 1)
    return true
end

function collect_perm!(v::AbstractVector)
    a = [copy(v)]
    while next_perm!(v) push!(a, copy(v)) end
    a
end

struct Perm{T<:AbstractVector} v::T end
perm(v::AbstractVector; kwargs...) = Perm(sort(v; kwargs...))

function Base.length(p::Perm)
    v = p.v
    N = factorial(length(v))
    i = firstindex(v)
    for j in eachindex(v)
        if isless(v[i], v[j])
            N = N ÷ factorial(j-i)
            i = j
        end
    end
    N = N ÷ factorial(lastindex(v)-i+1)
    N
end

Base.eltype(p::Perm{T}) where T<:AbstractVector = T

function Base.iterate(p::Perm)
    s = copy(p.v)
    t = copy(s)
    next_perm!(t) || return (s, nothing)
    u = copy(t)
    next_perm!(u) || return (s, (t, nothing))
    return (s, (t, u))
end

function Base.iterate(p::Perm, state)
    isnothing(state) && return nothing
    s, t = state
    isnothing(t) && return (s, nothing)
    u = copy(t)
    next_perm!(u) || return (s, (t, nothing))
    return (s, (t, u))
end

#function Base.iterate(p::Perm, s = copy(p.v))
#    isnothing(s) && return nothing
#    next_perm!(s) || return (s, nothing)
#    return (copy(s), s)
#end

end

# %%
using BenchmarkTools

# %%
p = O.perm([NaN, NaN, 2, 1])
_, state = iterate(p)

# %%
iterate(p, state)

# %%
p = O.perm([NaN, NaN, 2, 1])
@btime collect($p)

# %%
v = sort([NaN, NaN, 2, 1])
@btime O.collect_perm!($v)

# %%
w = sort(v)
@show w
while O.next_perm!(w)
    @show w
end

# %%
test_next_perm!(v) = (c = 1; while O.next_perm!(v) c += 1 end; c)
@show s = [1, 1, 2, 2, 2, 3, 4, 5, 6, 6]
@btime test_next_perm!($s)

# %%
c = @btime collect($(O.perm(s)))
length(c)

# %%
@code_warntype O.next_perm!(s)

# %%
using OffsetArrays
@show t = OffsetArray([1, 1, 2, 2, 2, 3, 4, 5, 6, 6], -5:4)
@btime collect($(O.perm(t)))

# %%
@btime O.collect_perm!($(sort(t)))

# %%
@code_warntype O.next_perm!(t)

# %%
