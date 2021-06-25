# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using Distributions, SpecialFunctions, QuadGK, Plots
pyplot(fmt = :svg)

α = range(0.1, 10; length=101)
# θ = 2

EX(α) = 2α
ElogX(α) = log(2) + digamma(α) # digamma(x) = ψ(x) = d(log Γ(x))/dx

F(α) = quadgk(x -> x*pdf(Gamma(α, 2.0), x), 0, Inf)[1]
G(α) = quadgk(x -> log(x)*pdf(Gamma(α, 2.0), x), 0, Inf)[1]

y = EX.(α)
z = ElogX.(α)
y_numint = F.(α)
z_numint = G.(α)

P = plot(α, y; label="E[X|α]", xlabel="α", legend=:topleft)
plot!(α, y_numint; label="num. int.", ls=:dash)

Q = plot(α, z; label="E[log X|α]", xlabel="α", legend=:bottomright)
plot!(α, z_numint; label="num. int.", ls=:dash)

plot(P, Q; size=(800, 300))

# %%
