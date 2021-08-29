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
using ConcreteStructs

@concrete struct Foo a; b; c; d; e; f; g; h; i; j; k; l; m; n; o; p; q; r; s; t; u; v; w; x; y; z end

names = fieldnames(Foo)
kwargs = Expr(:parameters, (Expr(:kw, name, v) for (v, name) in enumerate(names))...)
@eval Foo($kwargs) = Foo($(names...))

# %%
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

let foo = Foo(o = 9999)
    @unpackall_Foo foo
    @show a b c d e f g h i j k l m n o p q r s t u v w x y z
end;

# %%
using Parameters

let foo = Foo(o = 9999)
    @unpack a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z = foo
    @show a b c d e f g h i j k l m n o p q r s t u v w x y z
end;

# %%
methods(Foo)

# %%
?@unpackall_Foo

# %%
?@unpack

# %%
