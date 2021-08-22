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
Iterator generating tuples of repeat length of `v` and value `v` in `A`
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

# %%
