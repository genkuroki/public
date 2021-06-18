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

# %%
using Tricks

# %%
?static_hasmethod

# %%
struct A end
struct B end
abstract type FooBar end
struct Foo <: FooBar x::Int end
foomethod(foo::Foo) = nothing
struct Bar <: FooBar x::Int end

# %%
foo = Foo(23)
hasmethod(foomethod, Tuple{Foo})

# %%
attr_dynamic(fb::T) where T <: FooBar = hasmethod(foomethod, Tuple{T}) ? A() : B()
f_dynamic(fb::FooBar) = f(attr_dynamic(fb), fb)
f(::A, fb) = fb.x + 100
f(::B, fb) = fb.x - 100

# %%
f_dynamic(Foo(23))

# %%
f_dynamic(Bar(-23))

# %%
@code_typed f_dynamic(Foo(23))

# %%
foo = Foo(23)
static_hasmethod(foomethod, Tuple{Foo})

# %%
attr_static(fb::T) where T <: FooBar = static_hasmethod(foomethod, Tuple{T}) ? A() : B()
f_static(fb::FooBar) = f(attr_static(fb), fb)

# %%
f_static(Foo(23))

# %%
f_static(Bar(-23))

# %%
@code_typed f_static(Foo(23))

# %%
@code_typed f_static(Bar(-23))

# %%
foomethod(fb::Bar) = nothing

# %%
f_static(Bar(-23))

# %%
@code_typed f_static(Bar(-23))

# %%
