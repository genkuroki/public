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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
half(x::Int)::Int = x/2

# %%
half(2)

# %%
half(3)

# %%
@code_lowered half(2)

# %%
function f(n::Integer)::Vector{Real}
    randn(n)
end

# %%
@time sum(randn(10^8))

# %%
@time sum(f(10^8))

# %%
abstract type AbstractFoo{T} end

struct Foo{T} <: AbstractFoo{T} a::T end
double(x::Foo) = Foo(2x.a)

struct Bar{T} <: AbstractFoo{T} a::T end

# %%
double(Foo(1.23))

# %%
double(Bar(1.23))

# %%
double(x::AbstractFoo) = throw(MethodError(double, (x,)))

# %%
double(Bar(1.23))

# %%
multiply(x::Foo, y::Number) = Foo(x.a * y) 

# %%
multiply(Foo(3), 10.0)

# %%
multiply(Bar(3), 10.0)

# %%
multiply(x::AbstractFoo, y) = error("`myltiply` has not been implemented for $(typeof(x))")

# %%
multiply(Bar(3), 10.0)

# %%
multiply(Foo(3), [10.0])

# %%
multiply(x::AbstractFoo, y) = throw(MethodError(multiply, (x, y)))

# %%
multiply(Foo(3), [10.0])

# %%
function double!(a::Array)
    @. a = 2a
end

# %%
v = [1, 2, 3, 4]
@show double!(v)
@show v;

# %%
w = view(v, 2:3)

# %%
double!(w)

# %%
function double!(v)
    @. v = 2v
end

# %%
@show double!(w)
@show w
@show v;

# %%
