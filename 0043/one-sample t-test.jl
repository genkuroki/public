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
#     display_name: Julia 1.9.2
#     language: julia
#     name: julia-1.9
# ---

# %% tags=[]
using Distributions
using QuadGK
using Random
using StatsBase
using StatsPlots
default(fmt=:png, titlefontsize=9)

function sim(dist, n; μ=mean(dist), correction=true, L=10^6)
    D = correction ? TDist(n-1) : Normal()
    pval = Vector{Float64}(undef, L)
    Xtmp = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = Xtmp[tid]
        rand!(dist, X)
        X̄ = mean(X)
        S = std(X)
        T = √n * (X̄ - μ) / S
        pval[i] = 2ccdf(D, abs(T))
    end
    _ecdf_pval = ecdf(pval)
    ecdf_pval(x) = _ecdf_pval(x)
    ecdf_pval
end

function plot_sim(dist, n; μ=mean(dist), correction=true, L=10^6)
    testname = correction ? "t-test" : "no correction"
    distname = replace(string(dist), r"{[^}]*}"=>"")
    ecdf_pval = sim(dist, n; μ, correction, L)
    println("n = $n => coverage probability (α=5%) = ", 1 - ecdf_pval(0.05))
    plot(ecdf_pval, 0, 0.1; label="")
    plot!(identity; label="", c=:grey, ls=:dot)
    plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
    plot!(xguide="nominal significance level α", yguide="probability of P-value ≤ α")
    title!("$testname, $distname, n=$n")
    plot!(size=(400, 420))
end

mydistname(dist) = replace(string(dist), r"{[^}]*}"=>"")

function mydistname(dist::MixtureModel)
    d = mydistname.(dist.components)
    d = replace.(d, r", σ=1.0"=>"")
    p = string.(dist.prior.p)
    name = "$(p[1]) $(d[1])"
    name *= prod(" + $(p[i]) $(d[i])" for i in 2:length(d))
    name
end

function plot_sims(dist; ns=(5, 10, 20, 100), correction=true, L=10^6)
    distname = mydistname(dist)
    println("skewness of $distname = ", skewness(dist))
    println("kurtosis of $distname = ", kurtosis(dist))
    println()
    PP = []
    for n in ns
        P = plot_sim(dist, n; correction, L)
        push!(PP, P)
    end
    println()
    plot(PP...; size=(800, 840))
end

myskewness(dist) = skewness(dist)

function myskewness(dist::MixtureModel)
    μ, σ = mean(dist), std(dist)
    quadgk(x -> pdf(dist, x) * ((x - μ)/σ)^3, extrema(dist)...)[1]
end

mykurtosis(dist) = kurtosis(dist)

function mykurtosis(dist::MixtureModel)
    μ, σ = mean(dist), std(dist)
    quadgk(x -> pdf(dist, x) * ((x - μ)/σ)^4, extrema(dist)...)[1] - 3
end

function plot_sims2(dist; ns=[10; 20:20:200; 300:100:500], L=10^6, α=0.05, kwargs...)
    distname = mydistname(dist)
    println("skewness of $distname = ", myskewness(dist))
    println("kurtosis of $distname = ", mykurtosis(dist))
    println()
    xs = collect(ns)
    println("n ∈ $xs")
    println()
    
    μ, σ = mean(dist), std(dist)
    if σ == Inf
        μ, σ = 0.0, 1.2
    end
    P = plot(x -> pdf(dist, x), μ-5σ, μ+5σ; label="")
    title!(distname)
    plot!(; kwargs...)
    
    prob1 = similar(ns, Float64)
    prob2 = similar(prob1)
    for (i, n) in enumerate(ns)
        ecdf_pval1 = sim(dist, n; correction=true, L)
        prob1[i] = ecdf_pval1(α)
        ecdf_pval2 = sim(dist, n; correction=false, L)
        prob2[i] = ecdf_pval2(α)
    end
    Q = plot(xs, prob1; label="t-test", marker=:o)
    plot!(xs, prob2; label="no correction", ls=:dash, marker=:diamond)
    hline!([α]; label="", ls=:dot, c=:red)
    plot!(xguide="sample size n", yguide="probability of P-values ≤ $α")
    plot!(ytick=0:0.01:1)
    title!(distname)
    
    plot(P, Q; size=(640, 800), layout=(2,1))
end

# %%
plot_sims(Normal(); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(Normal(); ns=(5, 10, 20, 100))

# %%
plot_sims(Uniform(); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(Uniform(); ns=(5, 10, 20, 100))

# %%
plot_sims(Exponential(); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(Exponential(); ns=(5, 10, 20, 100))

# %%
plot_sims(Exponential(); ns=(100, 200, 500, 1000))

# %%
plot_sims(Beta(0.1, 0.1); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(Beta(0.1, 0.1); ns=(5, 10, 20, 100))

# %%
plot_sims(TDist(4); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(TDist(4); ns=(5, 10, 20, 100))

# %%
plot_sims(LogNormal(); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(LogNormal(); ns=(5, 10, 20, 100))

# %%
plot_sims(LogNormal(); ns=(100, 300, 1000, 3000), correction=false)

# %%
plot_sims(LogNormal(); ns=(100, 300, 1000, 3000))

# %%
plot_sims2(Normal())

# %%
plot_sims2(Uniform())

# %%
plot_sims2(Exponential(); ns=[20:20:200; 250; 500; 750; 1000])

# %%
plot_sims2(Beta(0.1, 0.1); ylim=(-0.2, 5))

# %%
plot_sims2(TDist(5))

# %%
plot_sims2(TDist(4))

# %%
plot_sims2(TDist(3))

# %%
plot_sims2(TDist(2))

# %%
plot_sims2(TDist(1.1))

# %%
plot_sims2(Gamma(2, 3))

# %%
mixnormal = MixtureModel([Normal(0, 1.5), Normal(2.5)], [0.7, 0.3])
plot_sims2(mixnormal)

# %%
mixnormal = MixtureModel([Normal(0, 1.5), Normal(5)], [0.7, 0.3])
plot_sims2(mixnormal)

# %%
mixnormal = MixtureModel([Normal(), Normal(10)], [0.9, 0.1])
plot_sims2(mixnormal; ns=[20:20:100; 150:50:250; 500; 750; 1000])

# %%
