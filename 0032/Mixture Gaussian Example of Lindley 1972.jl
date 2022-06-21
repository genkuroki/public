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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# Mixture Gaussian Example of Lindley 1972 in https://tminka.github.io/papers/minka-pathologies.pdf

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, titlefontsize=10)
using Optim

# %%
mixnormal(μ, σ) = MixtureModel([Normal(), Normal(μ, σ)], [0.5, 0.5])

# %%
Random.seed!(45105963)
n = 10
X = round.(rand(mixnormal(4, 1.5), n); digits=1)

# %%
f(μ, t) = loglikelihood(mixnormal(μ, exp(t)), X)

# %%
μ = range(-1, 8, 361)
t = range(-20, 2, 200)
z = f.(μ', t)
heatmap(μ, t, z)
title!("log likelihood")
plot!(; xguide="μ", yguide="log σ")

# %%
μ = range(-1, 8, 361)
t = range(-100, 2, 200)
z = f.(μ', t)
heatmap(μ, t, z)
title!("log likelihood")
plot!(; xguide="μ", yguide="log σ")

# %%
μ = range(-1, 9, 361)
t = range(-1, 2, 200)
z = exp.(f.(μ', t))
heatmap(μ, t, z)
title!("likelihood")
plot!(; xguide="μ", yguide="log σ")

# %%
@show o = optimize(x -> -f(x...), [4.0, 1.0])
o.minimizer[1], exp(o.minimizer[2])

# %%
for k in eachindex(X)
    o = optimize(x -> -f(x...), [X[k], -50])
    @show k, X[k], o.minimizer[1], exp(o.minimizer[2])
end

# %%
P = plot(; legend=:topleft)
for k in eachindex(X)
    plot!(t -> f(X[k], -t), -2, 70; label="$k", ls=:auto)
end
P

# %%
