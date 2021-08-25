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
#     display_name: Julia 1.7.0-beta3
#     language: julia
#     name: julia-1.7
# ---

# %%
#Base.findmin(f, X) = mapfoldl(x -> (f(x), x), min, X)

# %%
#ts = range(0, 2; length=201)
#findmin(sinpi, ts)

# %%
#rosenbrock2d((x, y),) = (1 - x)^2 + 100(y - x^2)^2
#xs = ys = range(-5, 5; length=1001)
#findmin(rosenbrock2d, Iterators.product(xs, ys))

# %%
struct Foo{F} f::F end
Base.findmin(foo::Foo, X) = mapfoldl(x -> (foo.f(x), x), min, X)

rosenbrock2d((x, y),) = (1 - x)^2 + 100(y - x^2)^2
xs = ys = range(-5, 5; length=1001)
findmin(Foo(rosenbrock2d), Iterators.product(xs, ys))

# %%
module My
findmin(f, X) = mapfoldl(x -> (f(x), x), min, X)
end

rosenbrock2d((x, y),) = (1 - x)^2 + 100(y - x^2)^2
xs = ys = range(-5, 5; length=1001)
My.findmin(rosenbrock2d, Iterators.product(xs, ys))

# %%
rosenbrock2d((x, y),) = (1 - x)^2 + 100(y - x^2)^2
xs = ys = range(-5, 5; length=1001)
findmin(rosenbrock2d, Iterators.product(xs, ys))

# %%
Base.keys(itr::Iterators.ProductIterator) = Iterators.product(keys.(itr.iterators)...)
findmin(rosenbrock2d, Iterators.product(xs, ys))

# %%
argmin(rosenbrock2d, Iterators.product(xs, ys))

# %%
valargmin(f, X) = (arg = argmin(f, X); (f(arg), arg))
valargmin(rosenbrock2d, Iterators.product(xs, ys))

# %%
Base.keys(itr::Iterators.ProductIterator) = Iterators.product(keys.(itr.iterators)...)
argmin(rosenbrock2d∘last, pairs(Iterators.product(xs, ys)))

# %%
Base.keys(itr::Iterators.ProductIterator) = Iterators.product(keys.(itr.iterators)...)
valargindmin(f, X) = ((ind, arg) = argmin(f∘last, pairs(X)); (f(arg), arg, ind))
valargindmin(rosenbrock2d, Iterators.product(xs, ys))

# %%
