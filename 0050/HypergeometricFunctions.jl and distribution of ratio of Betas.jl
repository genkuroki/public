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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# https://www.someweekendreading.blog/beta-ratios/

# %%
using Distributions
using HypergeometricFunctions
using QuadGK
using SpecialFunctions
using StatsFuns
using StatsPlots
default(fmt=:png)

function pdf_betaratio(κ, λ, μ, ν, x)
    if x < 0
        0.0
    elseif x ≤ 1
        exp(logbeta(κ+μ, ν) - logbeta(κ, λ) - logbeta(μ, ν) + xlogy(κ-1, x) + 
            log(max(0, _₂F₁(κ+μ, 1-λ, κ+μ+ν, x))))
    elseif x < Inf
        exp(logbeta(κ+μ, λ) - logbeta(κ, λ) - logbeta(μ, ν) + xlogy(-μ-1, x) +
        log(max(0, _₂F₁(κ+μ, 1-ν, κ+μ+λ, 1/x))))
    else
        Inf
    end
    
end

# %%
κ, λ, μ, ν = 2, 3, 4, 5
f = x -> pdf_betaratio(κ, λ, μ, ν, x)
@show quadgk(f, 0, Inf)[1]
L = 10^6
X = rand(Beta(κ, λ), L)
Y = rand(Beta(μ, ν), L)
R = @. X / Y
a, b = -0.2, min(5, maximum(R))
#a, b = 0.9, 1.1
plot()
stephist!(R; norm=true, label="")
plot!(f, extrema(R)...; label="")
plot!(xlim=(a, b))

# %%
κ, λ, μ, ν = 5, 10, 8, 13
f = x -> pdf_betaratio(κ, λ, μ, ν, x)
@show quadgk(f, 0, Inf)[1]
L = 10^6
X = rand(Beta(κ, λ), L)
Y = rand(Beta(μ, ν), L)
R = @. X / Y
a, b = -0.2, min(5, maximum(R))
#a, b = 0.9, 1.1
plot()
stephist!(R; norm=true, label="")
plot!(f, extrema(R)...; label="")
plot!(xlim=(a, b))

# %%
κ, λ, μ, ν = 20, 30, 15, 30
f = x -> pdf_betaratio(κ, λ, μ, ν, x)
@show quadgk(f, 0, Inf)[1]
L = 10^6
X = rand(Beta(κ, λ), L)
Y = rand(Beta(μ, ν), L)
R = @. X / Y
a, b = -0.2, min(5, maximum(R))
#a, b = 0.9, 1.1
plot()
stephist!(R; norm=true, label="")
plot!(f, extrema(R)...; label="")
plot!(xlim=(a, b))

# %%
κ, λ, μ, ν = 20, 30, 30, 30
f = x -> pdf_betaratio(κ, λ, μ, ν, x)
@show quadgk(f, 0, Inf)[1]
L = 10^6
X = rand(Beta(κ, λ), L)
Y = rand(Beta(μ, ν), L)
R = @. X / Y
a, b = -0.2, min(5, maximum(R))
#a, b = 0.9, 1.1
plot()
stephist!(R; norm=true, label="")
plot!(f, extrema(R)...; label="")
plot!(xlim=(a, b))

# %%
