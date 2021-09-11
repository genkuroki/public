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
sinpi(1/6)

# %%
apply_Function(f::Function, x) = f(x)

apply_Function(sinpi, 1/6)

# %%
struct Foo{F} f::F end
(foo::Foo)(x) = foo.f(x)

Foo(sinpi)(1/6)

# %%
apply_Function(Foo(sinpi), 1/6)

# %%
apply(f, x) = f(x)

apply(Foo(sinpi), 1/6)

# %%
