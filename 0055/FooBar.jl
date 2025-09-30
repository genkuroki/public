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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
module O

abstract type AbstractFoo end
plusone_a(foo::AbstractFoo) = foo.a += 1
double_b(foo::AbstractFoo) = foo.b *= 2

mutable struct Foo{T,U} <: AbstractFoo
    a::T
    b::U
end

abstract type AbstractBar <: AbstractFoo end
square_c(bar::AbstractBar) = bar.c ^= 2

mutable struct Bar{T,U,V} <: AbstractBar
    a::T
    b::U
    c::V
end

function Base.getproperty(foo::AbstractFoo, f::Symbol)
    f in fieldnames(typeof(foo)) && return getfield(foo, f)
    f_foo = Symbol(f, "_", string(foo))
    @eval function $f_foo(x...; y...) $f($foo, x...; y...) end
end

Base.propertynames(foo::T) where T<: AbstractFoo = [
    fieldnames(T)...;
    (m -> m.name).(methodswith(T; supertypes=true))
] |> unique

end

# %%
@show foo = O.Foo(1, 2.0)

@show O.plusone_a(foo)
@show foo

@show O.double_b(foo)
@show foo

@show bar = O.Bar(1, 2.0, "three")

@show O.plusone_a(bar)
@show bar

@show O.double_b(bar)
@show bar

@show O.square_c(bar)
@show bar
;

# %%
@show foo = O.Foo(1, 2.0)

@show foo.plusone_a()
@show foo

@show foo.double_b()
@show foo

@show bar = O.Bar(1, 2.0, "three")

@show bar.plusone_a()
@show bar

@show bar.double_b()
@show bar

@show bar.square_c()
@show bar
;

# %%
