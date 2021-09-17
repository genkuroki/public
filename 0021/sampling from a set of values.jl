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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
using StatsBase

sample([2, 3, 5], ProbabilityWeights([0.2, 0.3, 0.5]))

# %%
using StatsBase

a = sample([2, 3, 5], ProbabilityWeights([0.2, 0.3, 0.5]), 10^6)
[(x, count(==(x), a)/length(a)) for x in sort(unique(a))]

# %%
?sample

# %%
using Distributions

rand(DiscreteNonParametric([2, 3, 5], [0.2, 0.3, 0.5]))

# %%
using Distributions

a = rand(DiscreteNonParametric([2, 3, 5], [0.2, 0.3, 0.5]), 10^6)
[(x, count(==(x), a)/length(a)) for x in sort(unique(a))]

# %%
?DiscreteNonParametric

# %%
