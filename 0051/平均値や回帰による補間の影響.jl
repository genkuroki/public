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
#     display_name: Julia 1.11.1
#     language: julia
#     name: julia-1.11
# ---

# %%
using Distributions
using LinearAlgebra
using Random
Random.seed!(4649373)
using StatsPlots
default(fmt=:png, size=(400, 400), titlefontsize=10)

rd(x) = round(x; sigdigits=2)

function plot_ols(x, y; α=0.05, m=length(x),
        xlim=extrema(x), ms=3, ma=0.7, kwargs...)
    n = length(x)
    r = 2
    X = float(x) .^ (0:1)'
    β̂ = X \ y
    ŷ = X * β̂
    
    f̂(x) = evalpoly(x, β̂)
    ŝ = norm(y - ŷ)/√(n - r)
    xx(x) = [1, x]
    XX = X'X
    g(x) = ŝ * √dot(xx(x), XX \ xx(x))
    h(x) = ŝ * √(1 + dot(xx(x), XX \ xx(x)))
    t = quantile(TDist(n - r), 1 - α/2)
    
    plot(; legend=:topleft)
    scatter!(x[1:m], y[1:m]; label="data", c=1, ms, ma, msc=:auto)
    m != n && scatter!(x[m+1:n], y[m+1:n]; label="", c=5, ms, ma, msc=:auto)
    a, b = xlim
    a, b = a - 0.05(b-a), b + 0.05(b-a)
    plot!(f̂, a, b; label="", c=2, lw=1.3)
    plot!(x -> f̂(x) - t*g(x), a, b; label="$(100(1-α))% CI", c=3, ls=:dash)
    plot!(x -> f̂(x) + t*g(x), a, b; label="", c=3, ls=:dash)
    plot!(x -> f̂(x) - t*h(x), a, b; label="$(100(1-α))% PI", c=4, ls=:dashdot)
    plot!(x -> f̂(x) + t*h(x), a, b; label="", c=4, ls=:dashdot)
    title!("n=$n,  betahat=$(rd.(β̂)),  shat=$(rd(ŝ))")
    plot!(; kwargs...)
end

function plot_interpolations(;
        seed = 12345,
        β = [0, 1],
        σ = 0.3,
        n = 200,
        m = round(Int, 0.7n), 
        x = (Random.seed!(seed); 2 .+ rand(n)),
        u = σ * randn(n),
        y = evalpoly.(x, (β,)) + u,
        ylim = extrema(y) .* (0.95, 1.05),
        size = (800, 800),
        ms = 3,
    )
    @show β σ n m
    
    # オリジナルデータ
    P1 = plot_ols(x, y; xlim=extrema(x), ylim, ms)

    # ランダムな欠損
    P2 = plot_ols(x[1:m], y[1:m]; xlim=extrema(x), ylim, ms)

    # xが大きな所が欠損
    xsort = sort(1:n; by=i->x[i])
    P3 = plot_ols(x[xsort][1:m], y[xsort][1:m]; xlim=extrema(x), ylim, ms)

    # yが大きな所が欠損
    ysort = sort(1:n; by=i->y[i])
    P4 = plot_ols(x[ysort][1:m], y[ysort][1:m]; xlim=extrema(x), ylim, ms)
    
    # ランダムな欠損を平均値で補間
    xx = x[1:m]
    yy = y[1:m]
    xxx = x
    yyy = [yy; fill(mean(yy), n-m)]
    Q2 = plot_ols(xxx, yyy; xlim=extrema(x), ylim, m, ms)

    # xが大きな所の欠損を平均値で補間
    xsort = sort(1:n; by=i->x[i])
    xx = x[xsort][1:m]
    yy = y[xsort][1:m]
    xxx = x[xsort]
    yyy = [yy; fill(mean(yy), n-m)]
    Q3 = plot_ols(xxx, yyy; xlim=extrema(x), ylim, m, ms)

    # yが大きな所の欠損を平均値で補間
    ysort = sort(1:n; by=i->y[i])
    xx = x[ysort][1:m]
    yy = y[ysort][1:m]
    xxx = x[ysort]
    yyy = [yy; fill(mean(yy), n-m)]
    Q4 = plot_ols(xxx, yyy; xlim=extrema(x), ylim, m, ms)
    
    # ランダムな欠損を回帰で補間
    xx = x[1:m]
    yy = y[1:m]
    b = (float(xx) .^ (0:1)') \ yy
    xxx = x
    yyy = [yy; [evalpoly(x, b) for x in xxx[m+1:n]]]
    R2 = plot_ols(xxx, yyy; xlim=extrema(x), ylim, m, ms)
    
    # xが大きな所の欠損を回帰で補間
    xsort = sort(1:n; by=i->x[i])
    xx = x[xsort][1:m]
    yy = y[xsort][1:m]
    b = (float(xx) .^ (0:1)') \ yy
    xxx = x[xsort]
    yyy = [yy; [evalpoly(x, b) for x in xxx[m+1:n]]]
    R3 = plot_ols(xxx, yyy; xlim=extrema(x), ylim, m, ms)

    # yが大きな所の欠損を回帰で補間
    ysort = sort(1:n; by=i->y[i])
    xx = x[ysort][1:m]
    yy = y[ysort][1:m]
    b = (float(xx) .^ (0:1)') \ yy
    xxx = x[ysort]
    yyy = [yy; [evalpoly(x, b) for x in xxx[m+1:n]]]
    R4 = plot_ols(xxx, yyy; xlim=extrema(x), ylim, m, ms)
    
    plot(P1, P2, Q2, R2; size, layout=(2, 2)) |> display
    plot(P1, P3, Q3, R3; size, layout=(2, 2)) |> display
    plot(P1, P4, Q4, R4; size, layout=(2, 2)) |> display
end

# %%
plot_interpolations(; n=100)

# %%
plot_interpolations(; n=200)

# %%
plot_interpolations(; n=500, ms=2.5)

# %%
plot_interpolations(; n=1000, ms=2)

# %%
plot_interpolations(; seed=4649, n=300, m=240, x=rand(Normal(4, 0.5), 300), σ=1.0)

# %%
