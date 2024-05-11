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
using StatsPlots
default(fmt=:png)

# %%
?NegativeBinomial

# %% [markdown]
# * μ = mean(NegativeBinomial(r, p)) = (1-p)*r/p
# * σ² = var(NegativeBinomial(r, p)) = (1-p)*r/p^2
# * α = μ^2/σ² = (1-p)*r
# * θ = σ²/μ = 1/p

# %%
mypdf(disy, x) = pdf(dist, x)
mypdf(dist::DiscreteUnivariateDistribution, x) = pdf(dist, round(x))

function plot_with_gammadist(dist; 
        μ = mean(dist), σ² = var(dist), σ = var(dist),
        a = max(minimum(dist)-1, μ - 4σ), b = min(maximum(dist), μ + 4σ),
        kwargs...
    )
    @show dist
    @show μ
    @show σ²
    α, θ = μ^2/σ², σ²/μ
    @show gamma = Gamma(α, θ)
    plot(x -> mypdf(dist, x), a, b; label="dist")
    plot!(gamma, a, b; label="gamma", ls=:dash)
    plot!(; kwargs...)
end

# %%
plot_with_gammadist(NegativeBinomial(1, 0.1); xlim=(-2, 100), ylim=(-0.002, 0.12))

# %%
plot_with_gammadist(NegativeBinomial(1, 0.01); xlim=(-20, 1000), ylim=(-0.0002, 0.012))

# %%
plot_with_gammadist(NegativeBinomial(10, 0.6); xlim=(-0.5, 30))

# %%
plot_with_gammadist(NegativeBinomial(10, 0.2); xlim=(0, 150))

# %%
function plot_with_scaledgamma(dist ; 
        μ = mean(dist), σ² = var(dist),
        kwargs...
    )
    @show dist
    @show μ
    @show σ²
    α, θ = μ^2/σ², σ²/μ
    @show θ
    @show gamma = Gamma(α, 1)
    b = 7√α
    a = -1/θ
    plot(x -> mypdf(dist, θ*x)*θ, a, b; label="dist/θ")
    plot!(gamma, a, b; label="gamma", ls=:dash)
    plot!(; kwargs...)
end

# %%
plot_with_scaledgamma(NegativeBinomial(1, 0.1))

# %%
plot_with_scaledgamma(NegativeBinomial(1, 0.01))

# %%
plot_with_scaledgamma(NegativeBinomial(10, 0.6))

# %%
plot_with_scaledgamma(NegativeBinomial(10, 0.2))

# %%
