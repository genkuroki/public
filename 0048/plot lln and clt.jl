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
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, titlefontsize=12)

function plot_lln_clt(dist, n; L=10^6)
    @show dist
    @show n
    @show μ = mean(dist)
    @show σ = std(dist)
    X̄ = zeros(L)
    Xtmp = [zeros(eltype(dist), n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(dist, Xtmp[tid])
        X̄[i] = mean(X)
    end
    P = stephist(X̄; norm=true, label="sample means")
    vline!([μ]; label="", ls=:dash)
    plot!(xlim=(μ - 5σ, μ + 5σ))
    title!("law of large numbers")
    Q = stephist(X̄; norm=true, label="sample means")
    plot!(Normal(μ, σ/√n); label="", ls=:dash)
    plot!(xlim=(μ - 5σ/√n, μ + 5σ/√n))
    title!("central limit theorem")
    plot(P, Q; size=(600, 600), layout=(2, 1))
end

plot_lln_clt(Exponential(), 1000)

# %%
dist = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
plot(dist, -25, 25; label="dist", size=(600, 300))

# %%
plot_lln_clt(dist, 10)

# %%
plot_lln_clt(dist, 20)

# %%
plot_lln_clt(dist, 100)

# %%
plot_lln_clt(dist, 1000)

# %%
function plot_lln_clt_disc(dist, n=ntrials(dist); L=10^5)
    @show dist
    @show n
    @show μ = mean(dist)
    @show σ = std(dist)
    P = plot(x -> pdf(dist, round(x)),
        max(minimum(dist)-1, μ-5√n*σ), min(maximum(dist)+1, μ+5√n*σ); norm=true, label="")
    vline!([μ]; label="", ls=:dash)
    title!("law of large numbers")
    Q = plot(x -> pdf(dist, round(x)), μ - 5σ, μ + 5σ; norm=true, label="")
    plot!(Normal(μ, σ); label="", ls=:dash)
    title!("central limit theorem")
    plot(P, Q; size=(600, 600), layout=(2, 1))
end

plot_lln_clt_disc(Binomial(10000, 0.3))

# %%
plot_lln_clt_disc(Hypergeometric(30000, 70000, 10000), 10000)

# %%
plot_lln_clt_disc(BetaBinomial(10000, 3000, 7000))

# %%
plot_lln_clt_disc(Poisson(1000), 1000)

# %%
