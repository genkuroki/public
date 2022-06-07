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
using Random
using SpecialFunctions
using StatsPlots
default(fmt=:png, titlefontsize=10, size=(400, 250))

# %%
function unbiased_std(X)
    n = length(X)
    c = √((n-1)/2) * exp(loggamma((n-1)/2) - loggamma(n/2))
    c * std(X)
end

# %%
function plot_unbiased_std(; dist = Normal(4, 3), n = 10, L = 10^6, kwargs...)
    @show dist
    @show n
    @show σ = std(dist)
    tmp = zeros(n)
    unbiased_stds = [unbiased_std(rand!(dist, tmp)) for _ in 1:L]
    E_unbiased_stds = mean(unbiased_stds)
    @show E_unbiased_stds
    @show std(unbiased_stds)

    stephist(unbiased_stds; norm=true, label="\"unbiased\" std")
    vline!([E_unbiased_stds]; label="E[\"unbiased\" std]", c=:blue)
    vline!([σ]; label="std(dist)", c=:red, ls=:dash)
    plot!(; kwargs...)
end

# %%
plot_unbiased_std(dist=Normal(4, 3))

# %%
plot_unbiased_std(dist=Exponential(3), xtick=0:100, xlim=(-0.5, 10.5))

# %%
a = √(log((1+√(1+4*3^2))/2))
plot_unbiased_std(dist=LogNormal(0, a), xtick=0:100, xlim=(-0.5, 10.5))

# %%
plot_unbiased_std(dist=LogNormal(0, a), n=100, xtick=0:100, xlim=(-0.5, 10.5))

# %%
plot_unbiased_std(dist=LogNormal(0, a), n=1000, xtick=0:100, xlim=(-0.5, 10.5))

# %%
