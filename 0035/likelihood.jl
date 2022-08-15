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
using Random
Random.seed!(4649373)
using Plots

P(x, p) = p^x * (1 - p)^(1 - x)
L(X, p) = prod(P(x, p) for x in X)

@show X = rand(0:1, 12)
plot(p -> L(X, p), 0, 1; label="likelhood", xguide="x")

# %%
using Random
Random.seed!(4649373)
using Plots

sufficient_statistic(X) = (length(X), sum(X))

function likelihood(X, p)
    n, k = sufficient_statistic(X)
    p^k * (1 - p)^(n - k)
end

@show X = rand(0:1, 12)
plot(p -> likelihood(X, p), 0, 1; label="likelhood", xguide="x")

# %%
using BenchmarkTools

X = rand(0:1, 10^3)

a = @btime L($X, 0.5)
b = @btime likelihood($X, 0.5)
@show a b

# %%
