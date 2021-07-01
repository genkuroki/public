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
valargmax(f, X) = (x = argmax(f, X); (f(x), x))
valindargmax(f, X) = valargmax(f∘last, pairs(X))
valargmin(f, X) = (x = argmin(f, X); (f(x), x))
valindargmin(f, X) = valargmin(f∘last, pairs(X))

# %%
using Plots

# %%
domain = range(0, 2; length=2001)
m, (i, x) = valindargmax(sinpi, domain)
@show m, (i, x)
plot(domain, sinpi; label="sin(πx)")
scatter!([x], [m]; label="maximum")
plot!(; size=(400, 250))

# %%
m, (i, x) = valindargmin(sinpi, range(0, 2; length=2001))
@show m, (i, x)
plot(domain, sinpi; label="sin(πx)")
scatter!([x], [m]; label="minimum")
plot!(; size=(400, 250))

# %%
f(x, y) = (1 - x)^2 + 100(y - x^2)^2

X = Y = range(-5, 5; length=1001)
@show m, (x, y) = valargmin(Base.splat(f), Iterators.product(X, Y))

surface(X, Y, log∘f; color=:rainbow, clim = (-7.5, 10), camera=(30, 60))
scatter!([x], [y], [-7.5]; label="minimum", color=:cyan)
plot!(; xtick = -5:5, ytick = -5:5, zlim=(-7.5, 10))

# %%
heatmap(X, Y, log∘f; color=:rainbow, clim=(-7.5, 10))
scatter!([x], [y]; label="", color=:cyan)
plot!(; xlim=extrema(X), ylim=extrema(Y), size=(500, 400))
plot!(; xtick = -5:5, ytick = -5:5)

# %%
