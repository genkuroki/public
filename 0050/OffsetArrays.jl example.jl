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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %%
using OffsetArrays

module O
mutable struct Foo{T} a::T end
end

function f!(x)
    a = x.a
    for i in eachindex(a)
        a[i] *= 2
    end
    x
end

function iterf!(x, n)
    for _ in 1:n
        f!(x)
    end
    x
end

@show v = Float64.(1:4)
@show ov = OffsetVector(v, 0:3)
@show x = O.Foo(ov)

@time f!(x)
@time f!(x)
@time iterf!(x, 4)
@time iterf!(x, 4)

# %%
