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
using Distributions, StatsPlots, Roots

tstat(μ, X) = (mean(X) - μ)/√(var(X) / length(X))
pvalue(μ, X) = 2ccdf(TDist(length(X) - 1), abs(tstat(μ, X)))

function ci(X, α)
    X̄ = mean(X)
    c = quantile(TDist(length(X) - 1), 1 - α/2)
    d = √(var(X) / length(X))
    X̄ - c*d, X̄ + c*d
end

isininterval(x, int) = first(int) ≤ x ≤ last(int)
isinci(μ, X, α) = isininterval(μ, ci(X, α))
prob_trueisinci(dist, n, α; L=10^6) = mean(isinci(mean(dist), rand(dist, n), α) for _ in 1:L)

# %%
ci_roots(X, α) = find_zeros(μ -> log(pvalue(μ, X)) - log(α), -1e2, 1e2)

# %%
@show X = rand(10)
@show ci(X, 0.05)
@show ci_roots(X, 0.05);

# %%
isinci_roots(μ, X, α) = isininterval(μ, ci_roots(X, α))
prob_trueisinci_roots(dist, n, α; L=10^4) = mean(isinci_roots(mean(dist), rand(dist, n), α) for _ in 1:L)

# %%
@time prob_trueisinci_roots(Normal(), 20, 0.05)

# %%
