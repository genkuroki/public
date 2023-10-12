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
using LinearAlgebra
using Distributions
using StatsPlots
default(fmt=:png)
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

expeczero(dist) = dist - mean(dist)

standardized(dist) = (dist - mean(dist)) / std(dist)

# %%
x = range(-2, 2, 10)
β = [0, 1, 0, -0.2]
dist = Normal(0, 0.5)
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 10)
β = [0, 1, 0, -0.2]
dist = expeczero(Exponential(0.5))
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 20)
β = [0, 1, 0, -0.2]
dist = expeczero(Exponential(0.5))
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 40)
β = [0, 1, 0, -0.2]
dist = expeczero(Exponential(0.5))
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 10)
β = [0, 1, 0, -0.2]
dist = 0.5standardized(LogNormal())
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 20)
β = [0, 1, 0, -0.2]
dist = 0.5standardized(LogNormal())
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 40)
β = [0, 1, 0, -0.2]
dist = 0.5standardized(LogNormal())
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 80)
β = [0, 1, 0, -0.2]
dist = 0.5standardized(LogNormal())
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 160)
β = [0, 1, 0, -0.2]
dist = 0.5standardized(LogNormal())
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 320)
β = [0, 1, 0, -0.2]
dist = 0.5standardized(LogNormal())
d = 3
plot_sim(x, β, dist, d)

# %%
