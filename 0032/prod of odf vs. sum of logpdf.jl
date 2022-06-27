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
floatmax(Float64)

# %%
floatmin(Float64)

# %%
using Distributions
using StatsFuns
using StatsPlots
default(fmt=:png, size=(400, 250))

logisticmodel(a, b, x) = Bernoulli(logistic(a + b*x))
randlogistic(a, b, x) = @. rand(logisticmodel(a, b, x))

a, b, n = 1, 2, 3000
x = rand(Uniform(-3, 2), n)
y = randlogistic(a, b, x)
scatter(x, y; label="", ytick=0:1, msc=:auto, alpha=0.5, ms=2,
    size=(800, 100), ylim=(-0.3, 1.2))

# %%
prod(pdf(logisticmodel(a, b, x[i]), y[i]) for i in eachindex(x, y))

# %%
sum(logpdf(logisticmodel(a, b, x[i]), y[i]) for i in eachindex(x, y))

# %%
