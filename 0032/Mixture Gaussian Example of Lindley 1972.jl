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
n = 10
dist_true = mixnormal(4, 1.5)

# %%
plot(x -> pdf(dist_true, x), -5, 11; label="dist_true")

# %%
#Random.seed!(45105963)
#@show X = sort(round.(rand(dist_true, n); digits=1))
X = [-0.8, -0.1, 0.3, 0.8, 1.1, 1.9, 2.8, 3.8, 5.3, 7.4]
@show X;

# %%
loglik(μ, t) = loglikelihood(mixnormal(μ, exp(t)), X)
lik(μ, t) = exp(loglik(μ, t))

# %%
μ = range(-1, 8, 361)
t = range(-20, 2, 200)
heatmap(μ, t, loglik)
title!("log likelihood")
plot!(; xguide="μ", yguide="log σ")

# %%
μ = range(-1, 8, 361)
t = range(-100, 2, 200)
heatmap(μ, t, loglik)
title!("log likelihood")
plot!(; xguide="μ", yguide="log σ")

# %%
μ = range(-1, 9, 361)
t = range(-1, 2, 200)
heatmap(μ, t, lik; xlim=extrema(μ), ylim=extrema(t), colorbar=false)
title!("likelihood")
plot!(; xguide="μ", yguide="log σ")
scatter!([4], [log(1.5)]; c=:cyan, label="true param", legend=:bottomleft)

# %%
@show o = optimize(x -> -loglik(x...), [4.0, 1.0])
o.minimizer[1], exp(o.minimizer[2])

# %%
for k in eachindex(X)
    o = optimize(x -> -loglik(x...), [X[k], -50])
    @show k, X[k], o.minimizer[1], exp(o.minimizer[2])
end

# %%
P = plot(; legend=:topleft)
for k in eachindex(X)
    plot!(t -> loglik(X[k], -t), -2, 70; label="X[$k]", ls=:auto)
end
plot!(xguide="-t = -log σ", yguide="log likelihood")

# %%
@show X[1]
@eval @show lik($(X[1]), log(7.5e-28))
@show likelihood_of_true_param = lik(4, log(1.5))
@show lik(X[1], log(7.5e-28)) / likelihood_of_true_param;

# %%
