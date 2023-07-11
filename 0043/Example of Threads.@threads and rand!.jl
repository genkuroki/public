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
#     display_name: Julia 1.9.2
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, legend=false, xguide="x̄", yguide="s²",
    markersize=1, markerstrokecolor=:auto, markeralpha=0.3)

function randmeanvar!(dist, X)
    rand!(dist, X)
    mean(X), var(X)
end

function sim(dist::UnivariateDistribution, n; L=10^4)
    X̄ = Vector{Float64}(undef, L)
    S² = similar(X̄)
    X = Vector{Float64}(undef, n)
    for i in 1:L
        X̄[i], S²[i] = randmeanvar!(dist, X)
    end
    X̄, S²
end

function sim_threads(dist::UnivariateDistribution, n; L=10^4)
    X̄ = Vector{Float64}(undef, L)
    S² = similar(X̄)
    Xtmp = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = Xtmp[Threads.threadid()]
        X̄[i], S²[i] = randmeanvar!(dist, X)
    end
    X̄, S²
end

# %%
dist = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
n = 1000

# %%
@time X̄, S² = sim(dist, n)
@time X̄, S² = sim(dist, n)
@time X̄, S² = sim(dist, n)
scatter(X̄, S²)

# %%
@show Threads.nthreads()
@time X̄, S² = sim_threads(dist, n)
@time X̄, S² = sim_threads(dist, n)
@time X̄, S² = sim_threads(dist, n)
scatter(X̄, S²)

# %%
