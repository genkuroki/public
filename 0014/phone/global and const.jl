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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# https://discourse.julialang.org/t/how-to-correctly-define-and-use-global-variables-in-the-module-in-julia/65720/8

# %%
# The `::Float64`'s below are redundant.
module A
global a = 100.0::Float64
function f(x)::Float64
    a = x^2
    return a
end
println("f(6.0) = ", f(6.0), ",  a = ", a)
end;

# %%
# A is equivalent to the following.
module B
a = 100.0
function f(x)::Float64
    a = x^2
    return a
end
println("f(6.0) = ", f(6.0), ",  a = ", a)
end;

# %%
# Don't write return-type in general.
# https://github.com/JuliaLang/julia/blob/master/doc/src/manual/functions.md#return-type
module B′
a = 100.0
function f(x)
    a = x^2
    return a
end
println("f(6.0) = ", f(6.0), ",  a = ", a)
end;

# %%
# We need `global` for the immutable global variable `a` in for-loop.
module C
a = 100.0
function f(x)
    global a = x^2
    return a
end
println("f(6.0) = ", f(6.0), ",  a = ", a)
end;

# %%
# We don't need `global` for the content `a[]` of the mutable global variable `a`.
module D
a = Ref(100.0)
function f(x)
    a[] = x^2
    return a[]
end
println("f(6.0) = ", f(6.0), ",  a[] = ", a[])
end;

# %%
# D is type-unstable.  We need `const` for type stability.
module E
const a = Ref(100.0)
function f(x)
    a[] = x^2
    return a[]
end
println("f(6.0) = ", f(6.0), ",  a[] = ", a[])
end;

# %%
@code_warntype D.f(6.0)

# %%
@code_warntype E.f(6.0)

# %%
