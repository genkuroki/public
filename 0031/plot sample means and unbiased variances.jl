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
using Distributions, StatsPlots, Random, Base.Threads
default(fmt = :png)

function plot_samplemean_and_unbiasedvar(dist, n; L=10^4, kwargs...)
    distname = replace(string(dist), r"{[^}]*}"=>"")
    X̄ = zeros(L)
    S² = zeros(L)
    tmp = [zeros(eltype(dist), n) for _ in 1:nthreads()]
    @threads for i in eachindex(X̄, S²)
        X = rand!(dist, tmp[threadid()])
        X̄[i] = mean(X)
        S²[i] = var(X)
    end
    scatter(X̄, S²; label="", ms=2.5, ma=0.3, msw=0)
    title!("$distname, n=$n", kwargs...)
end

# %%
plot_samplemean_and_unbiasedvar(Normal(), 10)

# %%
plot_samplemean_and_unbiasedvar(Uniform(-1.732, 1.732), 10)

# %%
plot_samplemean_and_unbiasedvar(Exponential(), 10)

# %%
plot_samplemean_and_unbiasedvar(Exponential(), 100)

# %%
plot_samplemean_and_unbiasedvar(Exponential(), 1000)

# %%
plot_samplemean_and_unbiasedvar(Poisson(), 10)

# %%
plot_samplemean_and_unbiasedvar(Poisson(), 100)

# %%
plot_samplemean_and_unbiasedvar(Poisson(), 1000)

# %%
