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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using StatsBase
using Distributions
using Plots

X = randn(100)
plot(x -> ecdf(X)(x), -5, 5; label="ecdf", legend=:topleft)
plot!(x -> cdf(Normal(), x); label="Normal()")

# %%
