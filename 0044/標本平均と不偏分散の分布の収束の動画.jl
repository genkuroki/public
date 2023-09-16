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
using QuadGK
using Random
using StatsPlots
default(fmt=:png, size=(600, 400))

# %%
distname(dist) = replace(string(dist), r"{[^\}]*}"=>"")

function distname(dist::MixtureModel)
    d = distname.(dist.components)
    d = replace.(d, r", σ=1.0"=>"")
    p = string.(dist.prior.p)
    name = "$(p[1]) $(d[1])"
    name *= prod(" + $(p[i]) $(d[i])" for i in 2:length(d))
    name
end

function distname(dist::LocationScale)
    μ, σ, ρ = params(dist)
    m = μ != 0 ? "$μ + " : ""
    s = σ != 1 ? "$σ " : ""
    m * s * distname(ρ)
end

function stdmoment(dist::ContinuousUnivariateDistribution, k;
        μ = mean(dist),
        σ = std(dist),
        a = max(minimum(dist), μ - 20σ),
        b = min(maximum(dist), μ + 20σ)
    )
    f(x) = ((x - μ)/σ)^k * pdf(dist, x)
    quadgk(f, a, b)[1]
end

Skewness(dist) = skewness(dist)
Skewness(dist::MixtureModel; kwargs...) = stdmoment(dist, 3; kwargs... )

Kurtosis(dist) = kurtosis(dist)
Kurtosis(dist::MixtureModel; kwargs...) = stdmoment(dist, 4; kwargs... ) - 3

function samplemeanvar(dist, n; L=10^4)
    X̄ = Vector{Float64}(undef, L)
    S² = similar(X̄)
    Xtmp = [Vector{eltype(dist)}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist, Xtmp[Threads.threadid()])
        X̄[i] = mean(X)
        S²[i] = var(X)
    end
    X̄, S²
end

function plot_samplemeanvar(dist, n; L=10^4, kwargs...)
    X̄, S² = samplemeanvar(dist, n; L)
    scatter(X̄, S²; kwargs...)
end

function plot_samplemeanvar!(dist, n; L=10^4, kwargs...)
    X̄, S² = samplemeanvar(dist, n; L)
    scatter!(X̄, S²; kwargs...)
end

function gif_samplemeanvar(dist; xlim=:auto, ylim=:auto,
        ns = [round(Int, 10^s) for s in range(log(10, 10), log(10, 2000), 200)],
        kwargs...)
    @show distname(dist)
    @show μ = mean(dist)
    @show σ² = var(dist)
    @show Skewness(dist)
    @show Kurtosis(dist)
    @gif for n in [fill(ns[begin], 40); ns; fill(ns[end], 60)]
        plot_samplemeanvar(dist, n; label="", ms=2, msc=:auto, ma=0.2)
        scatter!([μ], [σ²]; label="", m=:star)
        plot!(; xlim, ylim)
        plot!(xguide="sample mean", yguide="unbiased sample variance")
        title!("$(distname(dist)), n = $n", titlefontsize=11)
        plot!(; kwargs...)
    end
end

# %%
dist1 = Gamma(2, 3)
dist2 = MixtureModel([dist1], [1.0])
@show skewness(dist1) Skewness(dist2)
@show kurtosis(dist1) Kurtosis(dist2);

# %%
dist3 = MixtureModel([Normal(), Normal(50)], [0.98, 0.02])
@show Skewness(dist3)
@show Kurtosis(dist3);

# %%
dist = Gamma(2, 3)
n = 100
@show μ = mean(dist)
@show σ² = var(dist)
@show Skewness(dist)
@show Kurtosis(dist)
plot_samplemeanvar(dist, n; label="", ms=2, msc=:auto, ma=0.2)
scatter!([μ], [σ²]; label="", m=:star)
plot!(xguide="sample mean", yguide="unbiased sample variance")
title!("$(distname(dist)), n = $n", titlefontsize=11)

# %%
gif_samplemeanvar(Normal(); xlim=(-0.7, 0.7), ylim=(0, 2))

# %%
gif_samplemeanvar(Uniform(-1.7320508, 1.7320508); xlim=(-0.7, 0.7), ylim=(0, 2))

# %%
gif_samplemeanvar(Laplace(0, 0.70710678); xlim=(-0.7, 0.7), ylim=(0, 2))

# %%
gif_samplemeanvar(TDist(4); xlim=(-1, 1), ylim=(0, 6))

# %%
gif_samplemeanvar(Exponential(); xlim=(0.3, 1.7), ylim=(0, 3))

# %%
gif_samplemeanvar(Gamma(2, 3); xlim=(3, 9), ylim=(0, 50))

# %%
gif_samplemeanvar(Poisson(); xlim=(0.3, 1.7), ylim=(0, 3))

# %%
dist = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
gif_samplemeanvar(dist; xlim=(-1, 6), ylim=(0, 100))

# %%
plot(x -> pdf(dist, x), -4, 24; label="")
plot!(xguide="x", yguide="probability density")
title!("$(distname(dist))")

# %%
dist = MixtureModel([Normal(), Normal(50)], [0.98, 0.02])
gif_samplemeanvar(dist; xlim=(-1, 6), ylim=(0, 250))

# %%
plot(x -> pdf(dist, x), -4, 54; label="")
plot!(xguide="x", yguide="probability density")
title!("$(distname(dist))")

# %%
dist = MixtureModel([Normal(), Normal(4)], [0.9, 0.1])
gif_samplemeanvar(dist; xlim=(-1, 2), ylim=(0, 7))

# %%
plot(x -> pdf(dist, x), -4, 8; label="")
plot!(xguide="x", yguide="probability density")
title!("$(distname(dist))")

# %%
