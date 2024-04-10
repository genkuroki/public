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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10)

distname(dist) = replace(string(dist), r"{[^}]*}"=>"", r".="=>"")
distname(dist::InverseGamma) = "InverseGamma$(params(dist))"
function distname(dist::MixtureModel)
    n = ncomponents(dist)
    c = components(dist)
    p = probs(dist)
    "$(p[1]) " * distname(c[1]) * prod(" + $(p[k]) " * distname(c[k]) for k in 2:n)
end

function plot_clt2(dist, n; L=10^4, kwargs...)
    X̄ = zeros(L)
    S² = zeros(L)
    for i in 1:L
        X = rand(dist, n)
        X̄[i] = mean(X)
        S²[i] = var(X) # unbiased variance
    end
    
    scatter(X̄, S²; label="", ms=2, msc=:auto, ma=0.2, kwargs...)
    plot!(xguide="sample mean", yguide="unbiased sample variance")
    title!(distname(dist) * ", n = $n")
end

# %%
plot_clt2(Normal(), 10)

# %%
plot_clt2(Normal(), 10; yscale=:log10)

# %%
plot_clt2(Normal(), 100; yscale=:log10)

# %%
plot_clt2(Uniform(), 10; yscale=:log10)

# %%
plot_clt2(Uniform(), 100; yscale=:log10)

# %%
dist = Exponential()
@show μ = mean(dist)
@show σ = std(dist)
@show skewness(dist)
@show kurtosis(dist)
for k in 2:10
    n = 2^k
    plot_clt2(dist, n; yscale=:log10) |> display
end

# %%
dist = InverseGamma(4.01, 1)
@show μ = mean(dist)
@show σ = std(dist)
@show skewness(dist)
@show kurtosis(dist)
for k in 2:10
    n = 2^k
    plot_clt2(dist, n; yscale=:log10) |> display
end

# %%
dist = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
@show μ = mean(dist)
@show σ = std(dist)
#@show skewness(dist)
#@show kurtosis(dist)
for k in 2:10
    n = 2^k
    plot_clt2(dist, n) |> display
end

# %%
dist = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
@show μ = mean(dist)
@show σ = std(dist)
#@show skewness(dist)
#@show kurtosis(dist)
for k in 2:10
    n = 2^k
    plot_clt2(dist, n; yscale=:log10) |> display
end

# %%
dist = Poisson(1)
@show μ = mean(dist)
@show σ = std(dist)
@show skewness(dist)
@show kurtosis(dist)
for k in 2:10
    n = 2^k
    plot_clt2(dist, n) |> display
end

# %%
