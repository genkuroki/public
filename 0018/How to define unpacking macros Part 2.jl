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

# %% [markdown]
# * https://discourse.julialang.org/t/for-monte-carlo-simulation-with-same-code-same-algorithm-how-fast-is-julia-compared-with-fortran/67021/41
# * https://discourse.julialang.org/t/for-monte-carlo-simulation-with-same-code-same-algorithm-how-fast-is-julia-compared-with-fortran/67021/43
# * https://github.com/genkuroki/public/blob/main/0018/%40defkwargs.ipynb
# * https://github.com/genkuroki/public/blob/main/0018/%40defunpack.ipynb
# * https://github.com/genkuroki/public/blob/main/0018/How%20to%20define%20unpacking%20macros.ipynb
# * https://github.com/genkuroki/public/blob/main/0018/How%20to%20define%20unpacking%20macros%20Part%202.ipynb

# %%
using ConcreteStructs

@concrete struct Foo a; b; c; d; e; f; g; h; i; j; k; l; m; n; o; p; q; r; s; t; u; v; w; x; y; z end

names = fieldnames(Foo)
kwargs = Expr(:parameters, (Expr(:kw, name, v) for (v, name) in enumerate(names))...)
@eval Foo($kwargs) = Foo($(names...))

foo = Foo(a = "meow", m = 'm', z = 99.99)

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

let
    @unpackall_Foo foo
    @show a b c d e f g h i j k l m n o p q r s t u v w x y z
end;

# %%
using Parameters

let
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
