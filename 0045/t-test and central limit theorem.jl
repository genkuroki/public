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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

# %%
function plot_sample_means(;
        distx = Exponential(),
        disty = Exponential(),
        m = 10,
        n = 10, 
        L = 10^6
    )
    @show distx disty m n

    @show Δμ = mean(distx) - mean(disty)
    @show se = √(var(distx)/m + var(disty)/n)

    X̄ = zeros(L)
    Ȳ = zeros(L)

    Threads.@threads for i in 1:L
        X = rand(distx, m)
        Y = rand(disty, n)
        X̄[i] = mean(X)
        Ȳ[i] = mean(Y)
    end

    P1 = plot()
    stephist!(X̄; norm=true, label="size-$m sample mean of distx")
    plot!(Normal(mean(distx), std(distx)/√m); label="normal approx.", ls=:dash)

    P2 = plot()
    stephist!(Ȳ; norm=true, label="size-$n sample mean of disty")
    plot!(Normal(mean(disty), std(disty)/√n); label="normal aporox.", ls=:dash)

    P3 = plot()
    stephist!(X̄ - Ȳ; norm=true, label="diff. of sample means")
    plot!(Normal(Δμ, se); label="normal aporox.", ls=:dash)

    plot(P1, P2, P3; size=(600, 700), layout=(3, 1))
end

# %%
plot_sample_means(distx=Exponential(1), disty=Exponential(1), m=10, n=10)

# %%
plot_sample_means(distx=Exponential(1), disty=Exponential(1), m=7, n=13)

# %%
plot_sample_means(distx=Exponential(2), disty=Exponential(1), m=10, n=10)

# %%
plot_sample_means(distx=Exponential(2), disty=Exponential(1), m=7, n=13)

# %%
plot_sample_means(distx=Exponential(2), disty=Exponential(1), m=13, n=7)

# %%
using LinearAlgebra
using Distributions
using StatsPlots
default(fmt=:png)
using StaticArrays
using Random

function sim(x, β, dist, d; L=10^6)
    f(x) = evalpoly(x, β)
    n, r = length(x), d+1
    y0 = f.(x)
    A = x .^ (0:d)'
    ginvA = (A'A)\A'
    utmp = [zeros(n) for _ in 1:Threads.nthreads()]
    y = [zeros(n) for _ in 1:Threads.nthreads()]
    b = [zeros(r) for _ in 1:Threads.nthreads()]
    β̂ = zeros(r, L)
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        u = rand!(dist, utmp[tid])
        @. y[tid] = y0 + u
        mul!(b[tid], ginvA, y[tid])
        β̂[:, i] .= b[tid]
    end
    β̂
end

rd(x; sigdigits=4) = round(x; sigdigits)

function plot_sim(x, β, dist, d; L=10^6)
    f(x) = evalpoly(x, β)
    n, r = length(x), d+1
    y0 = f.(x)
    A = x .^ (0:d)'
    
    β̂ = sim(x, β, dist, d; L)
    
    PP = []
    
    @show dist mean(dist) std(dist)
    println()
    P = plot(dist; label="residual dist.")
    push!(PP, P)
    
    u = rand(dist, n)
    y = y0 + u
    P = scatter(x, y; label="example of data", msc=:auto, ms=3)
    plot!(f, extrema(x)...; label="f(x)", ls=:dash)
    push!(PP, P)
    
    @show β
    @show meanβ̂ = rd.(vec(mean(β̂; dims=2)))
    @show rd.(std(dist) * .√diag(inv(A'A)))
    @show stdβ̂ = rd.(vec(std(β̂; dims=2)))
    println()
    
    for i in 1:r
        P = stephist(β̂[i, :]; norm=true, label="betahat[$(i)]")
        plot!(Normal(meanβ̂[i], stdβ̂[i]); label="normal approx.", ls=:dash)
        push!(PP, P)
    end

    @show n d
    
    h = 1 + (r+1) ÷ 2
    plot(PP...; size=(800, 250h), layout=(h, 2))
end

# %%
x = π * (-1:0.1:1)
β = [0, 1, 0, -0.2]
_dist = Gamma(2, 0.5)
dist = _dist - mean(_dist)
d = 3

plot_sim(x, β, dist, d; L=10^6)

# %%
x = π * (-1:0.1:1)
β = [0, 1, 0, -0.2]
_dist = Gamma(2, 0.5)
dist = _dist - mean(_dist)
d = 5

plot_sim(x, β, dist, d; L=10^6)

# %%
x = π * (-1:0.1:1)
β = [0, 1, 0, -0.2]
_dist = Gamma(2, 0.5)
dist = _dist - mean(_dist)
d = 7

plot_sim(x, β, dist, d; L=10^6)

# %%
x = π * (-1:0.5:1)
β = [0, 1, 0, -0.2]
_dist = Gamma(2, 0.5)
dist = _dist - mean(_dist)
d = 3

plot_sim(x, β, dist, d; L=10^6)

# %%
x = π * (-1:0.2:1)
β = [0, 1, 0, -0.2]
_dist = Gamma(2, 0.5)
dist = _dist - mean(_dist)
d = 3

plot_sim(x, β, dist, d; L=10^6)

# %%
x = π * (-1:0.2:1)
β = [0, 1, 0, -0.2]
_dist = Exponential(0.5)
dist = _dist - mean(_dist)
d = 3

plot_sim(x, β, dist, d; L=10^6)

# %%
x = π * (-1:0.2:1)
β = [0, 1, 0, -0.2]
_dist = Normal(0, 0.5)
dist = _dist - mean(_dist)
d = 3

plot_sim(x, β, dist, d; L=10^6)

# %%
