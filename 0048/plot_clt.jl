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

function plot_clt(dist, n; L=10^5, kwargs...)
    X̄ = zeros(L)
    for i in 1:L
        X = rand(dist, n)
        X̄[i] = mean(X)
    end
    histogram(X̄; norm=true, alpha=0.3, label="sample means", kwargs...)
    plot!(Normal(mean(dist), std(dist)/√n); label="Normal(μ, σ/√n)")
    title!(distname(dist) * ", n = $n")
end

# %%
using Distributions
using Random
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

function plot_clt(dist, n; L=10^5, kwargs...)
    X̄ = zeros(L)
    Xtmp = [zeros(n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(dist, Xtmp[tid])
        X̄[i] = mean(X)
    end
    histogram(X̄; norm=true, alpha=0.3, label="sample means", kwargs...)
    plot!(Normal(mean(dist), std(dist)/√n); label="Normal(μ, σ/√n)")
    title!(distname(dist) * ", n = $n")
end

# %%
plot_clt(Normal(1, 2), 10)

# %%
plot_clt(Bernoulli(0.3), 10; bin=(0.5:10.5)/10)

# %%
plot_clt(Bernoulli(0.3), 20; bin=(0.5:20.5)/20)

# %%
plot_clt(Bernoulli(0.3), 30; bin=(0.5:30.5)/30)

# %%
plot_clt(Uniform(), 1)

# %%
plot_clt(Uniform(), 2)

# %%
plot_clt(Uniform(), 3)

# %%
plot_clt(Uniform(), 4)

# %%
plot_clt(Uniform(), 5)

# %%
plot_clt(Uniform(), 6)

# %%
dist = Exponential()
@show μ = mean(dist)
@show σ = std(dist)
@show skewness(dist)
@show kurtosis(dist)
for k in 0:10
    n = 2^k
    plot_clt(dist, n; xlim=(μ-5σ/√n, μ+5σ/√n)) |> display
end

# %%
dist = InverseGamma(4.01, 1)
@show μ = mean(dist)
@show σ = std(dist)
@show skewness(dist)
@show kurtosis(dist)
for k in 0:10
    n = 2^k
    plot_clt(dist, n; xlim=(μ-5σ/√n, μ+5σ/√n)) |> display
end

# %%
dist = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
@show μ = mean(dist)
@show σ = std(dist)
#@show skewness(dist)
#@show kurtosis(dist)
for k in 0:10
    n = 2^k
    plot_clt(dist, n; xlim=(μ-4σ/√n, μ+8σ/√n)) |> display
end

# %%
