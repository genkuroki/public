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
# But `collect` causes memory allocations.

# %%
@time collect(XtimesY);

# %%
axes(XtimesY)

# %%
CartesianIndices(axes(XtimesY))

# %%
Base.keys(pr::Iterators.ProductIterator) = CartesianIndices(axes(pr))
val, idx = findmax(XtimesY) do (x, y); f(x, y) end

# %%
keys(XtimesY)
@time keys(XtimesY);

# %%
X[idx[1]], Y[idx[2]]

# %%
pt = argmax(XtimesY) do (x, y); f(x, y) end

# %%
valargmax(f, X) = (x = argmax(f, X); (f(x), x))
valargmax(XtimesY) do (x, y); f(x, y) end

# %%
valargmax(X) = valargmax(Base.Fix1(getindex, X), keys(X))
valargmax((pt -> f(pt...)).(XtimesY))

# %%
argmax((pt -> f(pt...)).(XtimesY))

# %%
findmax(pairs(X)) do (i, x); sin(x) end

# %%
findmax(sin, pairs(X))

# %%
pairs(pairs(X)) == pairs(X)

# %%
@which findmax(sin, pairs(X))

# %%
@which pairs(pairs(X))

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/b570546b68de16bd208ca76a20c1385919de18d6/base/iterators.jl#L236
#
# ```julia
# # pairs(v::Pairs) = v # listed for reference, but already defined from being an AbstractDict
# ```

# %%
?Base.Pairs

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/87e08d94a99ea17e7f72493947595eb6b63e0f09/base/essentials.jl#L32

# %%
Base.keys(rv::Iterators.Reverse) = keys(reverse(rv.itr))
findmax(sin, Iterators.Reverse(X))

# %%
Base.keys(en::Iterators.Enumerate) = keys(en.itr)
findmax(Iterators.enumerate(X)) do (i, x); sin(x) end

# %%
Base.keys(zp::Iterators.Zip) = Base.OneTo(length(zp))
findmax(Iterators.zip(X, reverse(Y))) do (x, y); cos(x)*sin(y) end

# %%
Base.keys(ac::Iterators.Accumulate) = keys(ac.itr)
findmax(sin, Iterators.accumulate(+, X; init=π))

# %%
Base.keys(tk::Iterators.Take) = Base.OneTo(tk.n)
findmax(sin, Iterators.take(reverse(X), 100))

# %%
methods(keys, Main)

# %%
?keys

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/16f433bb13cfc87eea21d26a797dac7b34a41d86/base/abstractdict.jl#L71

# %%
