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

@concrete struct Foo a; b; c end

macro with_kwargs(typename::Symbol)
    T = Core.eval(__module__, typename)
    names = fieldnames(T)
    kwargs = Expr(:parameters, names...)
    :($typename($kwargs) = $typename($(names...))) |> esc
end

macro with_kwargs(typename::Symbol, kwargs_expr)
    T = Core.eval(__module__, typename)
    names = fieldnames(T)
    D = Core.eval(__module__, kwargs_expr)
    kwargs = Expr(:parameters)
    for name in names
        val = get(D, name, nothing)
        if isnothing(val)
            push!(kwargs.args, name)
        else
            push!(kwargs.args, Expr(:kw, name, val))
        end
    end
    :($typename($kwargs) = $typename($(names...))) |> esc
end


# %%
@macroexpand @with_kwargs Foo

# %%
@with_kwargs Foo

# %%
methods(Foo)

# %%
Foo(a = 1, b = 2.0, c = "three")

# %%
@macroexpand @with_kwargs Foo (a = 1, c = "three")

# %%
@with_kwargs Foo (a = 1, c = "three")

# %%
methods(Foo)

# %%
Foo(b = 2.0)

# %%
default = (a = 1, b = 2.0, c = "three")
@macroexpand @with_kwargs Foo default

# %%
@with_kwargs Foo default

# %%
methods(Foo)

# %%
Foo()

# %%
