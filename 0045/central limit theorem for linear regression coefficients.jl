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
    
    PP = []
    
    @show dist mean(dist) std(dist)
    println()
    P = plot(dist; label="residual dist.")
    push!(PP, P)
    
    u = rand(dist, n)
    y = y0 + u
    b = A \ y
    f_lse(x) = evalpoly(x, b)
    P = scatter(x, y; label="example of data", msc=:auto, ms=3)
    plot!(f, extrema(x)...; label="true")
    plot!(f_lse, extrema(x)...; label="LSE", ls=:dash)
    push!(PP, P)
    
    β̂ = sim(x, β, dist, d; L)

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
ECDF(A, x) = count(≤(x), A)/length(A)

function sim_pval(x, β, dist, d; β₀ = β, L=10^6)
    f(x) = evalpoly(x, β)
    n, r = length(x), d+1
    y0 = f.(x)
    A = x .^ (0:d)'
    ginvA = (A'A)\A'
    sqrtdiaginvAA = .√diag(inv(A'A))
    nth = Threads.nthreads()
    utmp = [zeros(n) for _ in 1:nth]
    y = [zeros(n) for _ in 1:nth]
    b = [zeros(r) for _ in 1:nth]
    β̂ = zeros(r, L)
    uhat = [zeros(n) for _ in 1:nth]
    SEhat = zeros(r, L)
    T = zeros(r, L)
    pval = zeros(r, L)
    #Threads.@threads 
    for i in 1:L
        tid = Threads.threadid()
        u = rand!(dist, utmp[tid])
        @. y[tid] = y0 + u
        mul!(b[tid], ginvA, y[tid])
        @. β̂[:, i] = b[tid]
        @. uhat[tid] = y[tid] - evalpoly(x, (b[tid],))
        shat = √(dot(uhat[tid], uhat[tid])/(n - r))
        @. SEhat[:, i] = shat * sqrtdiaginvAA
        @. T[:, i] = (b[tid] - β₀) / @view(SEhat[:, i])
        @. pval[:, i] = 2ccdf(TDist(n-r), abs(@view(T[:, i])))
    end
    β̂, SEhat, T, pval
end

function plot_sim_pval(x, β, dist, d; β₀ = β, α = 0.05, L=10^6)
    f(x) = evalpoly(x, β)
    n, r = length(x), d+1
    y0 = f.(x)
    A = x .^ (0:d)'
    
    PP = []
    
    @show dist mean(dist) std(dist)
    println()
    P = plot(dist; label="residual dist.")
    push!(PP, P)
    
    u = rand(dist, n)
    y = y0 + u
    b = A \ y
    f_lse(x) = evalpoly(x, b)
    P = scatter(x, y; label="example of data", msc=:auto, ms=3)
    plot!(f, extrema(x)...; label="true")
    plot!(f_lse, extrema(x)...; label="LSE", ls=:dash)
    push!(PP, P)
     
    β̂, SEhat, T, pval = sim_pval(x, β, dist, d; β₀, L)

    @show β
    @show meanβ̂ = rd.(vec(mean(β̂; dims=2)))
    @show rd.(std(dist) * .√diag(inv(A'A)))
    @show stdβ̂ = rd.(vec(std(β̂; dims=2)))
    println()
    
    αs = exp.(range(log(0.002), log(1), 400))
    _tick = Any[0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
    xtick = ytick = (_tick, string.(_tick))
    
    @show n d
    @show β β₀
    @show α
    for i in 1:r
        @show ECDF(pval[i, :], α)
        P = plot(αs, α -> ECDF(pval[i, :], α); label="betahat[$i]")
        plot!(αs, identity; label="", ls=:dash)
        plot!(αs, x -> 1.2x; label="", ls=:dot, c=:black, alpha=0.4)
        plot!(αs, x -> 0.8x; label="", ls=:dot, c=:black, alpha=0.4)
        plot!(; xscale=:log10, yscale=:log10, xtick, ytick, xrotation=90)
        plot!(xguide="α", yguide="P(P-value ≤ α)")
        push!(PP, P)
    end
    
    h = 1 + (r+1) ÷ 2
    plot(PP...; size=(640, 320h), layout=(h, 2))
    plot!(leftmargin=4Plots.mm)
    plot!(tickfontsize=6)
end

# %%
x = range(-2, 2, 10)
β = [0, 1, 0, -0.2]
dist = Normal(0, 0.5)
d = 3
plot_sim(x, β, dist, d)

# %%
x = range(-2, 2, 10)
β = [0, 1, 0, -0.2]
dist = Normal(0, 0.5)
d = 5
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
x = range(-2, 2, 10)
β = [0, 1, 0, -0.2]
dist = Normal(0, 0.5)
d = 3

plot_sim_pval(x, β, dist, d; L=10^6)

# %%
x = range(-2, 2, 10)
β = [0, 1, 0, -0.2]
dist = Normal(0, 0.5)
d = 3

plot_sim_pval(x, β, dist, d; L=10^6, β₀=zero(β))

# %%
x = range(-2, 2, 10)
β = [0, 1, 0, -0.2]
dist = expeczero(Exponential(0.5))
d = 3

plot_sim_pval(x, β, dist, d; L=10^6)

# %%
x = range(-2, 2, 10)
β = [0, 1, 0, -0.2]
dist = expeczero(Exponential(0.5))
d = 3

plot_sim_pval(x, β, dist, d; L=10^6, β₀=zero(β))

# %%
x = range(-2, 2, 20)
β = [0, 1, 0, -0.2]
dist = expeczero(Exponential(0.5))
d = 3

plot_sim_pval(x, β, dist, d; L=10^6)

# %%
x = range(-2, 2, 20)
β = [0, 1, 0, -0.2]
dist = expeczero(Exponential(0.5))
d = 3

plot_sim_pval(x, β, dist, d; L=10^6, β₀=zero(β))

# %%
x = range(-2, 2, 40)
β = [0, 1, 0, -0.2]
dist = dist = expeczero(Exponential(0.5))
d = 3

plot_sim_pval(x, β, dist, d; L=10^6)

# %%
x = range(-2, 2, 40)
β = [0, 1, 0, -0.2]
dist = dist = expeczero(Exponential(0.5))
d = 3

plot_sim_pval(x, β, dist, d; L=10^6, β₀=zero(β))

# %%
x = range(-2, 2, 80)
β = [0, 1, 0, -0.2]
dist = dist = expeczero(Exponential(0.5))
d = 3

plot_sim_pval(x, β, dist, d; L=10^6)

# %%
x = range(-2, 2, 80)
β = [0, 1, 0, -0.2]
dist = dist = expeczero(Exponential(0.5))
d = 3

plot_sim_pval(x, β, dist, d; L=10^6, β₀=zero(β))

# %%
x = range(-2, 2, 40)
β = [0, 1, 0, -0.2]
dist = dist = 0.5standardized(LogNormal())
d = 3

plot_sim_pval(x, β, dist, d; L=10^6)

# %%
x = range(-2, 2, 80)
β = [0, 1, 0, -0.2]
dist = dist = 0.5standardized(LogNormal())
d = 3

plot_sim_pval(x, β, dist, d; L=10^6)

# %%
x = range(-2, 2, 160)
β = [0, 1, 0, -0.2]
dist = dist = 0.5standardized(LogNormal())
d = 3

plot_sim_pval(x, β, dist, d; L=10^6)

# %%
x = range(-2, 2, 320)
β = [0, 1, 0, -0.2]
dist = dist = 0.5standardized(LogNormal())
d = 3

plot_sim_pval(x, β, dist, d; L=10^6)

# %%
