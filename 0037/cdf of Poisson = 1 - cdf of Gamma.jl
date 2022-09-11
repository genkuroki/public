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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %%
using Distributions

λ, α, τ = 3, 5, 5
@show cdf(Poisson(λ*τ), α-1)
@show ccdf(Gamma(α, 1/τ), λ);

# %%
using Distributions
using StatsPlots

α, τ = 5, 2
plot(λ -> cdf(Poisson(λ*τ), α-1), 0, 8;
    label="cdf(Poisson(λτ), α-1)")
plot!(λ -> ccdf(Gamma(α, 1/τ), λ), 0, 8; ls=:dash,
    label="ccdf(Gamma(α, 1/τ), λ)")
title!("α = $α,  τ = $τ")

# %%
