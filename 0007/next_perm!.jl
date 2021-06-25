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

"""
zero allocation version of AquaIndigo's code
"""
function next_perm_old!(v::AbstractVector)
    length(v) ≤ 1 && return false    
    k = findlast(isless(v[i], v[i+1]) for i in firstindex(v):lastindex(v)-1)
    isnothing(k) && (reverse!(v); return false)
    i = k + firstindex(v) - 1    
    j = findlast(isless(v[i], v[j]) for j in i+1:lastindex(v)) + i
    v[i], v[j] = v[j], v[i]
    reverse!(v, i + 1)
    return true
end

"""
    next_perm!(v::AbstractVector)

changes `v` into the next permutation of `v` under the lexicographic order and returns `false` if it is the last permutation and `true` otherwise.

This is a variant of `Combinatorics.nextpermutation(m, t, state)`.

Examples

```
julia> v = [1, 2, 3]; println(v); while O.next_perm!(v) println(v) end
[1, 2, 3]
[1, 3, 2]
[2, 1, 3]
[2, 3, 1]
[3, 1, 2]
[3, 2, 1]
```

```
julia> v = [1, 2, 2]; println(v); while O.next_perm!(v) println(v) end
[1, 2, 2]
[2, 1, 2]
[2, 2, 1]
```
"""
function next_perm!(v::AbstractVector)
    length(v) ≤ 1 && return false
    i = lastindex(v) - 1
    @inbounds while i ≥ firstindex(v) && !isless(v[i], v[i+1]) i -= 1 end
    i < firstindex(v) && (reverse!(v); return false)
    j = lastindex(v)
    @inbounds while j > i && !isless(v[i], v[j]) j -= 1 end
    @inbounds v[i], v[j] = v[j], v[i]
    reverse!(v, i + 1)
    return true
end

function collect_perm!(v::AbstractVector)
    a = [copy(v)]
    while next_perm!(v) push!(a, copy(v)) end
    a
end

"""
    Perm(v::AbstractVector)

constructs the lexicographic order iterator of the all permutations of `v`.

Examples

```
julia> for s in O.Perm([1, 2, 3]) println(s) end
[1, 2, 3]
[1, 3, 2]
[2, 1, 3]
[2, 3, 1]
[3, 1, 2]
[3, 2, 1]
```

```
julia> for s in O.Perm([1, 2, 2]) println(s) end
[1, 2, 2]
[2, 1, 2]
[2, 2, 1]
```
"""
struct Perm{T<:AbstractVector}
    v::T
    Perm(v::T) where T<: AbstractVector = new{T}(sort(v))
end

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
@doc O.next_perm!

# %%
@doc O.Perm

# %%
v = [1, 2, 3]; println(v); while O.next_perm!(v) println(v) end

# %%
v = [1, 2, 2]; println(v); while O.next_perm!(v) println(v) end

# %%
for s in O.Perm([1, 2, 3]) println(s) end

# %%
for s in O.Perm([1, 2, 2]) println(s) end

# %%
v = sort([NaN, NaN, 2, 1])
@show v
while O.next_perm!(v)
    @show v
end

# %%
p = O.Perm([NaN, NaN, 2, 1])
for s in p
    @show s
end

# %%
using BenchmarkTools

# %%
p = O.Perm([NaN, NaN, 2, 1])
_, state = iterate(p)

# %%
iterate(p, state)

# %%
p = O.Perm([NaN, NaN, 2, 1])
@btime collect($p)

# %%
v = sort([NaN, NaN, 2, 1])
@btime O.collect_perm!($v)

# %%
test_next_perm!(v) = (c = 1; while O.next_perm!(v) c += 1 end; c)
@show s = [1, 1, 2, 2, 2, 3, 4, 5, 6, 6]
@btime test_next_perm!($s)

# %%
test_next_perm_old!(v) = (c = 1; while O.next_perm_old!(v) c += 1 end; c)
@show s = [1, 1, 2, 2, 2, 3, 4, 5, 6, 6]
@btime test_next_perm_old!($s)

# %%
c = @btime collect($(O.Perm(s)))
length(c)

# %%
@code_warntype O.next_perm!(s)

# %%
using OffsetArrays
@show t = OffsetArray([1, 1, 2, 2, 2, 3, 4, 5, 6, 6], -5:4)
@btime collect($(O.Perm(t)))

# %%
@btime O.collect_perm!($(sort(t)))

# %%
@code_warntype O.next_perm!(t)

# %%
