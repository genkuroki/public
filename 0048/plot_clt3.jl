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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %%
using LinearAlgebra
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10)

distname(dist) = replace(string(dist), r"{[^}]*}"=>"", r".="=>"")
distname(dist::InverseGamma) = "InverseGamma$(params(dist))"
function distname(dist::MixtureModel)
    n = ncomponents(dist)
    c = components(dist)
    p = probs(dist)
    "$(p[1]) " * distname(c[1]) * prod(" + $(p[k]) " * distname(c[k]) for k in 2:n)
end

function plot_clt3(dist, a=1:10; L=10^5, kwargs...)
    @show a
    n = length(a)
    X̄ = zeros(L)
    for i in 1:L
        X = rand(dist, n)
        X̄[i] = dot(a, X)
    end
    histogram(X̄; norm=true, alpha=0.3, label="dot(a, X)", kwargs...)
    plot!(Normal(sum(a)*mean(dist), norm(a)*std(dist)); label="normal approx.")
    title!(distname(dist) * ", n = $n")
end

# %%
plot_clt3(Uniform(), 1:10)

# %%
plot_clt3(Uniform(), (1:10).^2)

# %%
plot_clt3(Exponential(), 1:10)

# %%
plot_clt3(Exponential(), 1:100)

# %%
plot_clt3(Exponential(), round.(rand(100); digits=2))

# %%
plot_clt3(Exponential(), round.(randn(100); digits=2))

# %%
plot_clt3(MixtureModel([Normal(), Normal(20)], [0.95, 0.05]), 1:20)

# %%
plot_clt3(MixtureModel([Normal(), Normal(20)], [0.95, 0.05]), 1:100)

# %%
plot_clt3(MixtureModel([Normal(), Normal(20)], [0.95, 0.05]), 1:200)

# %%
plot_clt3(MixtureModel([Normal(), Normal(20)], [0.95, 0.05]), 1:1000)

# %%
plot_clt3(MixtureModel([Normal(), Normal(20)], [0.95, 0.05]), round.(rand(100); digits=2))

# %% tags=[]
plot_clt3(MixtureModel([Normal(), Normal(20)], [0.95, 0.05]), round.(randn(100); digits=2))

# %%
