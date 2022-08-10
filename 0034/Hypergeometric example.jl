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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png, size=(400, 250), titlefontsize=10)

k, N, n = 7, 10, 2
dist = Hypergeometric(k, N-k, n)
@show mean(dist), var(dist)
@show rationalize.((mean(dist), var(dist)))

xs = support(dist)
bar(xs, x -> pdf(dist, x); label="", title="$dist", xtick=xs)

# %%
using Distributions
using StatsPlots
default(fmt=:png, size=(400, 250), titlefontsize=10)

k, N, n = 20, 100, 10
dist = Hypergeometric(k, N-k, n)
@show mean(dist), var(dist)
@show rationalize.((mean(dist), var(dist)))

xs = support(dist)
bar(xs, x -> pdf(dist, x); label="", title="$dist", xtick=xs)

# %%
