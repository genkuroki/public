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

# %%
VERSION

# %%
f(x, y) = (1 - x)^2 + 100(y - x^2)^2
X = range(-5, 5, 501)[1:end-1]
Y = range(0, 5, 501)

# %%
findmin(sin, Iterators.reverse(collect(X)))

# %%
findmin(sin∘last, enumerate(X))

# %%
findmin(sin∘first, zip(X, Y))

# %%
findmin(sin, Iterators.accumulate(+, X))

# %%
findmin(sin, Iterators.take(X, 486))

# %%
findmin(v -> f(v...), Iterators.product(X, Y))

# %%
findmin(sin, Iterators.flatten((X, Y)))

# %%
Base.keys(rv::Iterators.Reverse) = keys(reverse(rv.itr))
Base.keys(en::Iterators.Enumerate) = keys(en.itr)
Base.keys(::Iterators.Zip) = Iterators.countfrom()
Base.keys(::Iterators.Filter) = Iterators.countfrom()
Base.keys(ac::Iterators.Accumulate) = keys(ac.itr)
Base.keys(::Iterators.Rest) = Iterators.countfrom()
Base.keys(::Iterators.Count) = Iterators.countfrom()
Base.keys(::Iterators.Take) = Iterators.countfrom()
Base.keys(::Iterators.Drop) = Iterators.countfrom()
Base.keys(::Iterators.TakeWhile) = Iterators.countfrom()
Base.keys(::Iterators.DropWhile) = Iterators.countfrom()
Base.keys(::Iterators.Cycle) = Iterators.countfrom()
Base.keys(::Iterators.Repeated) = Iterators.countfrom()
Base.keys(pr::Iterators.ProductIterator) = CartesianIndices(axes(pr))
Base.keys(::Iterators.Flatten) = Iterators.countfrom()
Base.keys(::Iterators.PartitionIterator) = Iterators.countfrom()

# %%
findmin(sin, Iterators.reverse(collect(X)))

# %%
findmin(sin∘last, enumerate(X))

# %%
findmin(sin∘first, zip(X, Y))

# %%
findmin(sin, Iterators.accumulate(+, X))

# %%
findmin(sin, Iterators.take(X, 486))

# %%
findmin(v -> f(v...), Iterators.product(X, Y))

# %%
findmin(sin, Iterators.flatten((X, Y)))

# %%
