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
using LinearAlgebra
using ForwardDiff
using Symbolics
using StaticArrays
using BenchmarkTools

# %%
function f(x)
    @show x
    @show zero(x)
    @assert x > zero(x)
    sin(x)/x
end

# %%
ForwardDiff.derivative(f, 1.0)

# %%
ForwardDiff.Dual{ForwardDiff.Tag{typeof(f), Float64}}(1.0,1.0) > 0

# %%
module My

struct Foo{T}
    a::T
    b::T
end
    
end

# %%
methods(My.Foo)

# %%
methods(My.Foo{Int})

# %%
My.Foo(a::T, b::T) where T<:Real = My.Foo{T}(a, b)

# %%
methods(My.Foo)

# %%
@which My.Foo(1im, 2im)

# %%
@which My.Foo(1, 2)

# %%
G(x, A) = exp(-dot(x, A, x)/2)

# %%
struct Bar{F, G}
    f::F
    ∇f::G
end

# %%
function Bar(f)
    ∇f(x, param) = ForwardDiff.gradient(x -> f(x, param), x)
    Bar(f, ∇f)
end

# %%
bar = Bar(G)

# %%
A = SA[
    2 -1
    -1 2
]

# %%
bar.f(SA[1, 2], A)

# %%
bar.∇f(SA[1, 2], A)

# %%
@variables a[axes(A)...], x[axes(SA[1, 2])...]

# %%
aa = collect(a)

# %%
collect(aa) - aa

# %%
xx = collect(x)

# %%
_G = G(xx, aa) # |> expand |> simplify

# %%
_∇G = Symbolics.gradient(_G, xx)

# %%
_∇G / _G

# %%
∇G = build_function(_∇G, xx, vec(aa); expression=Val(false))[1]

# %%
∇G(SA[1, 2], A)

# %%
function Bar(f, x, param)
    x_sym = collect(x)
    param_sym = collect(param)
    _f = simplify(expand(f(x_sym, param_sym)))
    _∇f = Symbolics.gradient(_f, x_sym)
    ff = build_function(_f, vec(x_sym), vec(param_sym); expression=Val(false))
    ∇f = build_function(_∇f, vec(x_sym), vec(param_sym); expression=Val(false))[1]
    Bar(ff, ∇f)
end

# %%
@variables x[1:2] a[1:2, 1:2]
car = Bar(G, x, a)

# %%
x = SA[1, 2]
A = SA[2 -1; -1 2]
car.f(x, A)

# %%
car.∇f(x, A)

# %%
@btime G($x, $A)
@btime $bar.f($x, $A)
@btime $car.f($x, $A)

# %%
@btime $∇G($x, $A)
@btime $bar.∇f($x, $A)
@btime $car.∇f($x, $A)

# %%
dG(x, A) = (A + A')*x/2*G(x, A)
#dG2(x, A) = A*x*G(x, A)
@btime $dG2($x, $A)

# %%
dG(xx, aa)

# %%
