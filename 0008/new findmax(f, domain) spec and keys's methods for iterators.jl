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
#     display_name: Julia 1.7.0-beta2
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# [old version of findmax(f, domain)](
# https://github.com/JuliaLang/julia/pull/35316/files#diff-97ad8b900c30b39398c32d78754c404f4a6df9a153d45284b2c11ab54837deb9R797)
#
# ```julia
# findmax(f, domain) = mapfoldl(x -> (f(x), x), _rf_findmax, domain)
#  _rf_findmax((fm, m), (fx, x)) = isless(fm, fx) ? (fx, x) : (fm, m)
# ```
#
# [current version](https://github.com/JuliaLang/julia/pull/41076/files#diff-97ad8b900c30b39398c32d78754c404f4a6df9a153d45284b2c11ab54837deb9R803)
#
# ```julia
# findmax(f, domain) = mapfoldl( ((k, v),) -> (f(v), k), _rf_findmax, pairs(domain) )
#  _rf_findmax((fm, im), (fx, ix)) = isless(fm, fx) ? (fx, ix) : (fm, im)
# ```
#
# Ref. [`findmax` and friends: confusing behaviour to be introduced in 1.7](https://discourse.julialang.org/t/findmax-and-friends-confusing-behaviour-to-be-introduced-in-1-7/61904)

# %%
using OffsetArrays

# %%
VERSION

# %%
f(x, y) = cos(x)*sin(y) + 0.1(x - y)

X = range(-2, 2; length=401)
Y = range(-2, 2; length=401)
XtimesY = Iterators.product(X, Y)

# %%
val, idx = findmax(XtimesY) do (x, y); f(x, y) end

# %% [markdown]
# `pairs` function used in `findmax(f, domain)` function requires `keys` method.

# %% tags=[]
methods(keys)

# %%
keys(collect(XtimesY))

# %%
val, idx = findmax(collect(XtimesY)) do (x, y); f(x, y) end

# %% [markdown]
# But `collect` causes memory allocations and I do not feel like it.

# %%
axes(XtimesY)

# %%
CartesianIndices(axes(XtimesY))

# %%
Base.keys(pr::Iterators.ProductIterator) = CartesianIndices(axes(pr))
F((x, y)) = f(x, y)

@show val, idx = findmax(F, XtimesY)
@show f(X[idx[1]], Y[idx[2]])
@show argmax(F.(XtimesY));

# %%
valargmax(f, X) = (x = argmax(f, X); (f(x), x))
valargmax(X) = valargmax(Base.Fix1(getindex, X), keys(X))
F((x, y)) = f(x, y)

@show val, idx = findmax(F, XtimesY)
@show val, (X[idx[1]], Y[idx[2]])
@show argmax(F, XtimesY)
@show valargmax(F, XtimesY)
@show valargmax(F.(XtimesY))
@show findmax(F, XtimesY);

# %%
Base.keys(rv::Iterators.Reverse) = keys(reverse(rv.itr))

@show val, idx = findmax(sin, Iterators.Reverse(X))
@show sin(reverse(X)[idx]);

# %%
Base.keys(en::Iterators.Enumerate) = keys(en.itr)
G((i, x)) = sin(x)

@show val, idx = findmax(G, enumerate(X))
@show sin(X[idx]);

# %%
Base.keys(zp::Iterators.Zip) = Base.OneTo(length(zp))
H((x, y)) = cos(x) * sin(y)

@show val, idx = findmax(H, zip(X, reverse(Y)))
@show cos(X[idx]) * sin(reverse(Y)[idx]);

# %%
Base.keys(ac::Iterators.Accumulate) = keys(ac.itr)
A = OffsetArray(range(-2, 2; length=401), -200:200)

@show val, idx = findmax(sin, Iterators.accumulate(+, A; init=π))
@show sin(sum(A[begin:idx]) + π);

# %%
Base.keys(tk::Iterators.Take) = Base.OneTo(length(tk))

@show val, idx = findmax(sin, Iterators.take(reverse(X), 100))
@show sin(reverse(X)[idx]);

# %%
methods(keys, Main)

# %%
?keys

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/16f433bb13cfc87eea21d26a797dac7b34a41d86/base/abstractdict.jl#L71

# %%
?Base.isgreater

# %%
?findmax

# %%
findmax([8, 0.1, -9, pi])

# %%
A = [1.0 2; 3 4]
findmax(A, dims=1)

# %%
methods(findmax)

# %%
findmax(-(-4:5).^2)

# %%
findmax(x -> -x^2, -4:5)

# %%
