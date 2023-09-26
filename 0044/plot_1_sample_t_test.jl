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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

countecdf(A, x) = count(≤(x), A)/length(A)

function pvalue_1_sample_t_test(X, μ)
    n = length(X)
    X̄ = mean(X)
    SEhat = std(X)/√n
    t = (X̄ - μ)/SEhat
    2ccdf(TDist(n-1), abs(t))
end

dist = Normal(2, 3)
μ₀ = mean(dist)

P = plot(dist; label="dist")
vline!([μ₀]; label="μ₀")
title!("population distribution")

n = 20
niters = 10^6
pval = zeros(niters)
Threads.@threads for i in 1:niters
    X = rand(dist, n)
    pval[i] = pvalue_1_sample_t_test(X, μ₀)
end

α = 0.05
println("(probability of P-value ≤ $α) = ", countecdf(pval, α))

Q = plot(α -> countecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot)
plot!(xguide="α", yguide="probability of P-value ≤ α")
title!("ECDF of P-values (n = $n)")

plot(P, Q; size=(720, 360))

# %%
function plot_1_sample_t_test(; dist=Normal(2, 3), μ₀=mean(dist), n=20, α=0.05, niters=10^6)
    P = plot(dist; label="dist")
    vline!([μ₀]; label="μ₀")
    title!("population distribution")

    pval = zeros(niters)
    Threads.@threads for i in 1:niters
        X = rand(dist, n)
        pval[i] = pvalue_1_sample_t_test(X, μ₀)
    end
    
    println("(probability of P-value ≤ $α) = ", countecdf(pval, α))
    
    Q = plot(α -> countecdf(pval, α), 0, 0.1; label="")
    plot!(identity; label="", ls=:dot)
    plot!(xguide="α", yguide="probability of P-value ≤ α")
    title!("ECDF of P-values (n = $n)")

    plot(P, Q; size=(720, 360))
end

plot_1_sample_t_test(dist = Normal(2, 3), n = 20)

# %%
plot_1_sample_t_test(; dist = Gamma(20, 3), n = 20)

# %%
plot_1_sample_t_test(; dist = Gamma(2, 3), n = 20)

# %%
plot_1_sample_t_test(; dist = Gamma(2, 3), n = 200)

# %%
plot_1_sample_t_test(; dist = Normal(0.2, 1), μ₀ = 0, n = 20)

# %%
plot_1_sample_t_test(; dist = Normal(0.2, 1), μ₀ = 0, n = 50)

# %%
plot_1_sample_t_test(; dist = Normal(0.2, 1), μ₀ = 0, n = 100)

# %%
plot_1_sample_t_test(; dist = Normal(0.2, 1), μ₀ = 0, n = 200)

# %%
@show dist = Gamma(20, 3)
@show μ = mean(dist)
@show σ = std(dist)
plot_1_sample_t_test(; dist, μ₀ = μ, n = 20)

# %%
@show dist = Gamma(20, 3)
@show μ = mean(dist)
@show σ = std(dist)
plot_1_sample_t_test(; dist, μ₀ = μ + 0.2σ, n = 20)

# %%
@show dist = Gamma(20, 3)
@show μ = mean(dist)
@show σ = std(dist)
plot_1_sample_t_test(; dist, μ₀ = μ + 0.2σ, n = 50)

# %%
@show dist = Gamma(20, 3)
@show μ = mean(dist)
@show σ = std(dist)
plot_1_sample_t_test(; dist, μ₀ = μ + 0.2σ, n = 100)

# %%
@show dist = Gamma(20, 3)
@show μ = mean(dist)
@show σ = std(dist)
plot_1_sample_t_test(; dist, μ₀ = μ + 0.2σ, n = 200)

# %%
@show dist = Gamma(2, 3)
@show μ = mean(dist)
@show σ = std(dist)
plot_1_sample_t_test(; dist, μ₀ = μ, n = 20)

# %%
@show dist = Gamma(2, 3)
@show μ = mean(dist)
@show σ = std(dist)
plot_1_sample_t_test(; dist, μ₀ = μ, n = 100)

# %%
@show dist = Gamma(2, 3)
@show μ = mean(dist)
@show σ = std(dist)
plot_1_sample_t_test(; dist, μ₀ = μ + 0.2σ, n = 100)

# %%
@show dist = Gamma(2, 3)
@show μ = mean(dist)
@show σ = std(dist)
plot_1_sample_t_test(; dist, μ₀ = μ + 0.2σ, n = 200)

# %%
using Random
using StatsBase: ecdf

function plot_1_sample_t_test_revised(; dist=Normal(2, 3), μ₀=mean(dist), n=20, α=0.05, niters=10^6)
    P = plot(dist; label="dist", title="population distribution")
    vline!([μ₀]; label="μ₀")

    nth = Threads.nthreads()
    Xtmp = [zeros(n) for _ in 1:nth]
    pval = zeros(niters)
    Threads.@threads for i in 1:niters
        tid = Threads.threadid()
        X = rand!(dist, Xtmp[tid])
        pval[i] = pvalue_1_sample_t_test(X, μ₀)
    end
    
    _pval = ecdf(pval)
    ecdf_pval(x) = _pval(x)
    println("(probability of P-value ≤ $α) = ", ecdf_pval(α))
    
    Q = plot(α -> ecdf_pval(α), 0, 0.1; label="", title="ECDF of P-values (n = $n)")
    plot!(identity; label="", ls=:dot)
    plot!(xguide="α", yguide="probability of P-value ≤ α")

    plot(P, Q; size=(720, 360))
end

plot_1_sample_t_test_revised()

# %%
@time plot_1_sample_t_test(; n=1000)
@time plot_1_sample_t_test_revised(; n=1000)

# %%
