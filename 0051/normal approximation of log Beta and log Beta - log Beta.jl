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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

α, β = 20, 180
beta1 = Beta(α, β)

L = 10^6
P = rand(beta1, L)
logP = @. log(P)
normal = Normal(log(α/(α+β)), √(β/α/(α+β+1)))

density(logP; label="density of log Beta($α, $β)")
plot!(normal; label="normal approximation", ls=:dash)
plot!(legend=:outertop, legendfontsize=12)

# %%
using Distributions
using StatsPlots
default(fmt=:png)

α, β = 20, 180
γ, δ = 10, 190
beta1 = Beta(α, β)
beta2 = Beta(γ, δ)

L = 10^6
P = rand(beta1, L)
Q = rand(beta2, L)
logP = @. log(P)
logQ = @. log(Q)
logRR = logP - logQ
normal = Normal(log(α/(α+β)) - log(γ/(γ+δ)), √(β/α/(α+β+1) + δ/γ/(γ+δ+1)))

density(logRR; label="density of log Beta($α, $β) − log Beta($γ, $δ)")
plot!(normal; label="normal approximation", ls=:dash)
plot!(legend=:outertop, legendfontsize=12)

# %%
