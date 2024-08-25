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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png)

function sim(dist, n; L=10^7)
    NCT = zeros(L)
    Xtmp = zeros(n)
    for i in 1:L
        X = rand!(dist, Xtmp)
        X̄, S² = mean(X), var(X)
        NCT[i] = X̄ / √(S²/n)
    end
    NCT
end

function plot_sim(dist, n; L=10^7, legend=:outertop, kwargs...)
    μ, σ = mean(dist), std(dist)
    NCT = sim(dist, n; L)
    nct = NoncentralT(n-1, √n*μ/σ)
    @show dist n μ σ nct
    @show mean(NCT), std(NCT)
    @show mean(nct), std(nct)
    @show √n*μ/σ
    stephist(NCT; norm=true, label="(sample mean) / √((unbiased variance)/n)", legend, kwargs...)
    plot!(nct; label="noncentral t-distribution", ls=:dash)
end

plot_sim(Normal(20, 10), 5; xlim=(0, 20))

# %%
plot_sim(Uniform(), 100)

# %%
plot_sim(Laplace(2, 3), 5; xlim=(-5, 10))

# %%
plot_sim(LogNormal(), 100)

# %%
L = 10^6
a, b, θ = 3, 7, 5
X = rand(Gamma(a, θ), L)
Y = rand(Gamma(b, θ), L)
P = @. X / (X + Y)
stephist(P; norm=true, label="X / (X + Y)")
plot!(Beta(a, b); label="Beta(a, b)", ls=:dash)
plot!(legendfontsize=12)

# %%
