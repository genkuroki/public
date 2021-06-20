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
# * https://twitter.com/kdwkshh/status/1406620780457119745
# * https://discourse.julialang.org/t/ann-jettest-jl-advanced-testing-toolset-for-julia/63229
# * https://github.com/aviatesk/JETTest.jl

# %%
VERSION

# %%
using JETTest
using Test

module O

struct Foo
    a::Vector
end

function Base.sum(foo::Foo)
    a = foo.a
    s = zero(eltype(a))
    for x in a
        s += x
    end
    s
end

end

foo = O.Foo(randn(10^3))
@show sum(foo)

println("\n", "="^80, "\n")

@report_dispatch sum(foo)

println("\n", "="^80, "\n")

@code_warntype sum(foo)

println("\n", "="^80, "\n")

@testset "check type-stabilities" begin
    @test_nodispatch sum(foo)
end

# %%
using JETTest
using Test

module O

struct Foo
    a::Vector{Real}
end

function Base.sum(foo::Foo)
    a = foo.a
    s = zero(eltype(a))
    for x in a
        s += x
    end
    s
end

end

foo = O.Foo(randn(10^3))
@show sum(foo)

println("\n", "="^80, "\n")

@report_dispatch sum(foo)

println("\n", "="^80, "\n")

@code_warntype sum(foo)

println("\n", "="^80, "\n")

@testset "check type-stabilities" begin
    @test_nodispatch sum(foo)
end

# %%
using JETTest
using Test

module O

struct Foo
    a::Array{Float64}
end

function Base.sum(foo::Foo)
    a = foo.a
    s = zero(eltype(a))
    for x in a
        s += x
    end
    s
end

end

foo = O.Foo(randn(10^3))
@show sum(foo)

println("\n", "="^80, "\n")

@report_dispatch sum(foo)

println("\n", "="^80, "\n")

@code_warntype sum(foo)

println("\n", "="^80, "\n")

@testset "check type-stabilities" begin
    @test_nodispatch sum(foo)
end

# %%
using JETTest
using Test

module O

struct Foo
    a::Array{Float64}
end

function Base.sum(foo::Foo)
    a = foo.a
    s = zero(eltype(a))
    for x in a
        s += x
    end
    s
end

end

foo = O.Foo(randn(10^3))
@show sum(foo)

println("\n", "="^80, "\n")

@report_dispatch frame_filter = (sv -> sv.mod === @__MODULE__) sum(foo)

println("\n", "="^80, "\n")

@code_warntype sum(foo)

println("\n", "="^80, "\n")

@testset "check type-stabilities" begin
    @test_nodispatch frame_filter = (sv -> sv.mod === @__MODULE__) sum(foo)
end

# %%
using JETTest
using Test

module O

struct Foo
    a::Vector{Float64}
end

function Base.sum(foo::Foo)
    a = foo.a
    s = zero(eltype(a))
    for x in a
        s += x
    end
    s
end

end

foo = O.Foo(randn(10^3))
@show sum(foo)

println("\n", "="^80, "\n")

@report_dispatch sum(foo)

println("\n", "="^80, "\n")

@code_warntype sum(foo)

println("\n", "="^80, "\n")

@testset "check type-stabilities" begin
    @test_nodispatch sum(foo)
end

# %%
