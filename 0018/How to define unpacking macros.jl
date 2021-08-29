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
struct Foo{A, B, C} a::A; b::B; c::C end

"""
`@unpackall_Foo(obj)` unpacks all fields of the object `obj` of type `Foo`.
"""
macro unpackall_Foo(obj)
    names = fieldnames(Foo)
    Expr(:(=),
        Expr(:tuple, names...),
        Expr(:tuple, (:($obj.$name) for name in names)...)
    ) |> esc
end

@doc @unpackall_Foo

# %%
@macroexpand @unpackall_Foo foo

# %%
@unpackall_Foo Foo(1, 2.0, "three")

# %%
a, b, c

# %%
foo = Foo(1, 2.0, "three")

function f(foo::Foo)
    @unpackall_Foo foo
    @show a b c
    return
end

f(foo)

# %%
