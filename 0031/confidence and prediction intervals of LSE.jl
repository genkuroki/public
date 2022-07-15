# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Distributions
using LinearAlgebra
using Random
Random.seed!(4649373)
using StatsPlots
default(fmt=:png, size=(400, 280), titlefontsize=10)

# %%
function plot_least_squares_estimate(x, y, f...; α = 0.05, kwargs...)
    n = length(x)
    r = length(f)
    X = [f(x) for x in x, f in f]
    β̂ = X \ y
    ŷ = X * β̂
    
    f̂(x) = sum(c * f(x) for (c, f) in zip(β̂, f))
    ŝ = norm(y - ŷ)/√(n - r)
    xx(x) = [f(x) for f in f]
    invXX = inv(X'X)
    g(x) = ŝ * √(xx(x)'invXX*xx(x))
    h(x) = ŝ * √(1 + xx(x)'invXX*xx(x))
    t = quantile(TDist(n - r), 1 - α/2)

    plot(; legend=:topleft)
    scatter!(x, y; label="data", c=1, msc=:auto)
    a, b = extrema(x)
    a, b = a - 0.05(b-a), b + 0.05(b-a)
    plot!(f̂, a, b; label="", c=:red, lw=2)
    plot!(x -> f̂(x) - t*g(x), a, b; label="$(100(1-α))% CI", c=3, ls=:dash)
    plot!(x -> f̂(x) + t*g(x), a, b; label="", c=3, ls=:dash)
    plot!(x -> f̂(x) - t*h(x), a, b; label="$(100(1-α))% PI", c=4, ls=:dashdot)
    plot!(x -> f̂(x) + t*h(x), a, b; label="", c=4, ls=:dashdot)
    plot!(; kwargs...)
end

# %%
n = 30
f₀(x) = sin(x)
x = range(-π, π, n)
y = f₀.(x) + 0.3randn(n)
scatter(x, y; label="data", legend=:topleft, msc=:auto)
plot!(f₀, extrema(x)...; label="y = f₀(x)", lw=2)

# %%
plot_least_squares_estimate(x, y, (x->x^k for k in 0:0)...)

# %%
plot_least_squares_estimate(x, y, (x->x^k for k in 0:1)...)

# %%
plot_least_squares_estimate(x, y, (x->x^k for k in 0:3)...)

# %%
plot_least_squares_estimate(x, y, (x->x^k for k in 0:5)...)

# %%
plot_least_squares_estimate(x, y, (x->x^k for k in 0:11)...; ylim=(-2.5, 2.5))

# %%
function plot_pvalue_functions(x, y, f...; α = 0.05, kwargs...)
    n = length(x)
    r = length(f)
    X = [f(x) for x in x, f in f]
    β̂ = X \ y
    ŷ = X * β̂
    
    f̂(x) = sum(c * f(x) for (c, f) in zip(β̂, f))
    s = norm(y - ŷ)/√(n - r)
    xx(x) = [f(x) for f in f]
    invXX = inv(X'X)
    g(x) = s * √(xx(x)'invXX*xx(x))
    h(x) = s * √(1 + xx(x)'invXX*xx(x))
    G(x, y) = 2ccdf(TDist(n-r), abs(y - f̂(x))/g(x))
    H(x, y) = 2ccdf(TDist(n-r), abs(y - f̂(x))/h(x))
    t = quantile(TDist(n - r), 1 - α/2)

    a, b = extrema(x)
    a, b = a - 0.05(b-a), b + 0.05(b-a)
    xs = range(a, b, 400)
    c, d = extrema(y)
    c, d = c - 0.1t*(d-c), d + 0.1t*(d-c)
    ys = range(c, d, 400)
    
    P = plot(; legend=:topleft, colorbar=false)
    heatmap!(xs, ys, sqrt∘G)
    scatter!(x, y; label="data", c=:cyan, msc=1)
    plot!(x -> f̂(x) - t*g(x), a, b; label="$(100(1-α)) CI", c=:pink, ls=:dash)
    plot!(x -> f̂(x) + t*g(x), a, b; label="", c=:pink, ls=:dash)
    plot!(; xlim=(a, b), ylim=(c, d))
    title!("P-value functon of CI")
    
    Q = plot(; legend=:topleft, colorbar=false)
    heatmap!(xs, ys, sqrt∘H)
    scatter!(x, y; label="data", c=:cyan, msc=1)
    plot!(x -> f̂(x) - t*h(x), a, b; label="$(100(1-α)) PI", c=:pink, ls=:dash)
    plot!(x -> f̂(x) + t*h(x), a, b; label="", c=:pink, ls=:dash)
    plot!(; xlim=(a, b), ylim=(c, d))
    title!("P-value function of PI")
    
    plot(P, Q; size=(800, 300), layout=(1, 2), kwargs...)
end

# %%
plot_pvalue_functions(x, y, (x->x^k for k in 0:0)...)

# %%
plot_pvalue_functions(x, y, (x->x^k for k in 0:1)...)

# %%
plot_pvalue_functions(x, y, (x->x^k for k in 0:3)...)

# %%
plot_pvalue_functions(x, y, (x->x^k for k in 0:5)...)

# %%
plot_pvalue_functions(x, y, (x->x^k for k in 0:11)...)

# %%
N = 50
X = randn(N)
Y = X + randn(N)
scatter(X, Y; label="data", legend=:topleft, msc=:auto)

# %%
plot_pvalue_functions(X, Y, (x->x^k for k in 0:1)...)

# %%
function plot_3d_pvalue_function(x, y, f...; α = 0.05, camera=(-45, 60), kwargs...)
    n = length(x)
    r = length(f)
    X = [f(x) for x in x, f in f]
    β̂ = X \ y
    ŷ = X * β̂
    
    f̂(x) = sum(c * f(x) for (c, f) in zip(β̂, f))
    s = norm(y - ŷ)/√(n - r)
    xx(x) = [f(x) for f in f]
    invXX = inv(X'X)
    g(x) = s * √(xx(x)'invXX*xx(x))
    G(x, y) = 2ccdf(TDist(n-r), abs(y - f̂(x))/g(x))
    t = quantile(TDist(n - r), 1 - α/2)

    a, b = extrema(x)
    a, b = a - 0.05(b-a), b + 0.05(b-a)
    xs = range(a, b, 400)
    c, d = extrema(y)
    c, d = c - 0.1t*(d-c), d + 0.1t*(d-c)
    ys = range(c, d, 400)
    
    P = plot(; colorbar=false)
    surface!(xs, ys, (x,y) -> G(x,y); camera, alpha=0.9)
    scatter3d!(x, y, zeros(length(x)); label="", c=:cyan, msc=:auto)
    plot!(; size=(600, 500))
end

# %%
pyplot()
plot_3d_pvalue_function(X, Y, (x->x^k for k in 0:1)...)

# %%
pyplot()
@gif for t in 0:5:359
    plot_3d_pvalue_function(X, Y, (x->x^k for k in 0:1)...; camera=(t, 60))
end

# %%
gr()

# %%
@gif for t in 0:5:359
    plot_3d_pvalue_function(X, Y, (x->x^k for k in 0:1)...; camera=(t, 60))
end

# %%
