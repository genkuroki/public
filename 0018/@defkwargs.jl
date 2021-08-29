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

@concrete struct Foo a; b; c end

macro defkwargs(typename::Symbol)
    T = Core.eval(__module__, typename)
    names = fieldnames(T)
    kwargs = Expr(:parameters, names...)
    :($typename($kwargs) = $typename($(names...))) |> esc
end

macro defkwargs(typename::Symbol, kwargs_expr)
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
@macroexpand @defkwargs Foo

# %%
@defkwargs Foo

# %%
methods(Foo)

# %%
Foo(a = 1, b = 2.0, c = "three")

# %%
@macroexpand @defkwargs Foo (a = 1, c = "three")

# %%
@defkwargs Foo (a = 1, c = "three")

# %%
methods(Foo)

# %%
Foo(b = 2.0)

# %%
default = (a = 1, b = 2.0, c = "three")
@macroexpand @defkwargs Foo default

# %%
@defkwargs Foo default

# %%
methods(Foo)

# %%
Foo()

# %%
